! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_kstat_1_8
    use scifort_descriptive_kstat_engine, only : kstat_all_assumed_rank, &
        kstat_kind_stat, kstat_kind_var, kstat_reduce_assumed_rank, &
        kstat_vector, kstatvar_vector
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface kstat
        module procedure kstat_all_rank1
        module procedure kstat_axis_rank1
        module procedure kstat_all_rank2
        module procedure kstat_axis_rank2
        module procedure kstat_all_rank3
        module procedure kstat_axis_rank3
        module procedure kstat_all_rank4
        module procedure kstat_axis_rank4
        module procedure kstat_all_rank5
        module procedure kstat_axis_rank5
        module procedure kstat_all_rank6
        module procedure kstat_axis_rank6
        module procedure kstat_all_rank7
        module procedure kstat_axis_rank7
        module procedure kstat_all_rank8
        module procedure kstat_axis_rank8
    end interface kstat

    interface kstat_all
        module procedure kstat_all_rank1
        module procedure kstat_all_rank2
        module procedure kstat_all_rank3
        module procedure kstat_all_rank4
        module procedure kstat_all_rank5
        module procedure kstat_all_rank6
        module procedure kstat_all_rank7
        module procedure kstat_all_rank8
    end interface kstat_all

    interface kstat_keepdims
        module procedure kstat_keepdims_rank1
        module procedure kstat_keepdims_rank2
        module procedure kstat_keepdims_rank3
        module procedure kstat_keepdims_rank4
        module procedure kstat_keepdims_rank5
        module procedure kstat_keepdims_rank6
        module procedure kstat_keepdims_rank7
        module procedure kstat_keepdims_rank8
    end interface kstat_keepdims

    interface kstatvar
        module procedure kstatvar_all_rank1
        module procedure kstatvar_axis_rank1
        module procedure kstatvar_all_rank2
        module procedure kstatvar_axis_rank2
        module procedure kstatvar_all_rank3
        module procedure kstatvar_axis_rank3
        module procedure kstatvar_all_rank4
        module procedure kstatvar_axis_rank4
        module procedure kstatvar_all_rank5
        module procedure kstatvar_axis_rank5
        module procedure kstatvar_all_rank6
        module procedure kstatvar_axis_rank6
        module procedure kstatvar_all_rank7
        module procedure kstatvar_axis_rank7
        module procedure kstatvar_all_rank8
        module procedure kstatvar_axis_rank8
    end interface kstatvar

    interface kstatvar_all
        module procedure kstatvar_all_rank1
        module procedure kstatvar_all_rank2
        module procedure kstatvar_all_rank3
        module procedure kstatvar_all_rank4
        module procedure kstatvar_all_rank5
        module procedure kstatvar_all_rank6
        module procedure kstatvar_all_rank7
        module procedure kstatvar_all_rank8
    end interface kstatvar_all

    interface kstatvar_keepdims
        module procedure kstatvar_keepdims_rank1
        module procedure kstatvar_keepdims_rank2
        module procedure kstatvar_keepdims_rank3
        module procedure kstatvar_keepdims_rank4
        module procedure kstatvar_keepdims_rank5
        module procedure kstatvar_keepdims_rank6
        module procedure kstatvar_keepdims_rank7
        module procedure kstatvar_keepdims_rank8
    end interface kstatvar_keepdims

    public :: kstat
    public :: kstat_all
    public :: kstat_keepdims
    public :: kstatvar
    public :: kstatvar_all
    public :: kstatvar_keepdims

contains

    pure function kstat_all_rank1(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank1

    pure function kstat_axis_rank1(x, n, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        if (axis /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        order_value = 2
        if (present(n)) order_value = n
        value = kstat_vector(x, order_value, nan_policy)
    end function kstat_axis_rank1

    pure function kstat_keepdims_rank1(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: order_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        order_value = 2
        if (present(n)) order_value = n
        allocate(values(1))
        values(1) = kstat_vector(x, order_value, nan_policy)
    end function kstat_keepdims_rank1

    pure function kstatvar_all_rank1(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank1

    pure function kstatvar_axis_rank1(x, n, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        if (axis /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        order_value = 2
        if (present(n)) order_value = n
        value = kstatvar_vector(x, order_value, nan_policy)
    end function kstatvar_axis_rank1

    pure function kstatvar_keepdims_rank1(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: order_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        order_value = 2
        if (present(n)) order_value = n
        allocate(values(1))
        values(1) = kstatvar_vector(x, order_value, nan_policy)
    end function kstatvar_keepdims_rank1

    pure function kstat_all_rank2(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank2

    pure function kstat_axis_rank2(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: order_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank2

    pure function kstat_keepdims_rank2(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: order_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank2

    pure function kstatvar_all_rank2(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank2

    pure function kstatvar_axis_rank2(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: order_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank2

    pure function kstatvar_keepdims_rank2(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: order_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank2

    pure function kstat_all_rank3(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank3

    pure function kstat_axis_rank3(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer :: order_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank3

    pure function kstat_keepdims_rank3(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: order_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank3

    pure function kstatvar_all_rank3(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank3

    pure function kstatvar_axis_rank3(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer :: order_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank3

    pure function kstatvar_keepdims_rank3(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: order_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank3

    pure function kstat_all_rank4(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank4

    pure function kstat_axis_rank4(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer :: order_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank4

    pure function kstat_keepdims_rank4(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: order_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank4

    pure function kstatvar_all_rank4(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank4

    pure function kstatvar_axis_rank4(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer :: order_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank4

    pure function kstatvar_keepdims_rank4(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: order_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank4

    pure function kstat_all_rank5(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank5

    pure function kstat_axis_rank5(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer :: order_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank5

    pure function kstat_keepdims_rank5(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: order_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank5

    pure function kstatvar_all_rank5(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank5

    pure function kstatvar_axis_rank5(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer :: order_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank5

    pure function kstatvar_keepdims_rank5(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: order_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank5

    pure function kstat_all_rank6(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank6

    pure function kstat_axis_rank6(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank6

    pure function kstat_keepdims_rank6(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: order_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank6

    pure function kstatvar_all_rank6(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank6

    pure function kstatvar_axis_rank6(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank6

    pure function kstatvar_keepdims_rank6(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: order_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank6

    pure function kstat_all_rank7(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank7

    pure function kstat_axis_rank7(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank7

    pure function kstat_keepdims_rank7(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: order_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank7

    pure function kstatvar_all_rank7(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank7

    pure function kstatvar_axis_rank7(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank7

    pure function kstatvar_keepdims_rank7(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: order_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank7

    pure function kstat_all_rank8(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_stat, nan_policy)
    end function kstat_all_rank8

    pure function kstat_axis_rank8(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstat_axis_rank8

    pure function kstat_keepdims_rank8(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: order_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_stat, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstat_keepdims_rank8

    pure function kstatvar_all_rank8(x, n, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        integer :: order_value

        order_value = 2
        if (present(n)) order_value = n
        value = kstat_all_assumed_rank(x, order_value, kstat_kind_var, nan_policy)
    end function kstatvar_all_rank8

    pure function kstatvar_axis_rank8(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: order_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kstatvar_axis_rank8

    pure function kstatvar_keepdims_rank8(x, n, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: n !! k-statistic variance order, default 2
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: order_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        order_value = 2
        if (present(n)) order_value = n
        call kstat_reduce_assumed_rank(x, axis_value, order_value, kstat_kind_var, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function kstatvar_keepdims_rank8

end module scifort_descriptive_kstat_1_8
