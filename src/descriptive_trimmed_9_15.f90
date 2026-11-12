! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_trimmed_9_15
    use scifort_descriptive_trimmed_engine, only : trimmed_kind_max, trimmed_kind_mean, &
        trimmed_kind_min, trimmed_kind_sem, trimmed_kind_std, trimmed_kind_var, &
        trimmed_reduce_assumed_rank, trimmed_stat_all_assumed_rank, trimmed_stat_vector
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface tmean
        module procedure tmean_all_rank9
        module procedure tmean_all_rank10
        module procedure tmean_all_rank11
        module procedure tmean_all_rank12
        module procedure tmean_all_rank13
        module procedure tmean_all_rank14
        module procedure tmean_all_rank15
        module procedure tmean_axis_rank9
        module procedure tmean_axis_rank10
        module procedure tmean_axis_rank11
        module procedure tmean_axis_rank12
        module procedure tmean_axis_rank13
        module procedure tmean_axis_rank14
        module procedure tmean_axis_rank15
    end interface tmean

    interface tmean_all
        module procedure tmean_all_rank9
        module procedure tmean_all_rank10
        module procedure tmean_all_rank11
        module procedure tmean_all_rank12
        module procedure tmean_all_rank13
        module procedure tmean_all_rank14
        module procedure tmean_all_rank15
    end interface tmean_all

    interface tmean_keepdims
        module procedure tmean_keepdims_rank9
        module procedure tmean_keepdims_rank10
        module procedure tmean_keepdims_rank11
        module procedure tmean_keepdims_rank12
        module procedure tmean_keepdims_rank13
        module procedure tmean_keepdims_rank14
        module procedure tmean_keepdims_rank15
    end interface tmean_keepdims

    interface tvar
        module procedure tvar_rank9
        module procedure tvar_rank10
        module procedure tvar_rank11
        module procedure tvar_rank12
        module procedure tvar_rank13
        module procedure tvar_rank14
        module procedure tvar_rank15
    end interface tvar

    interface tvar_all
        module procedure tvar_all_rank9
        module procedure tvar_all_rank10
        module procedure tvar_all_rank11
        module procedure tvar_all_rank12
        module procedure tvar_all_rank13
        module procedure tvar_all_rank14
        module procedure tvar_all_rank15
    end interface tvar_all

    interface tvar_keepdims
        module procedure tvar_keepdims_rank9
        module procedure tvar_keepdims_rank10
        module procedure tvar_keepdims_rank11
        module procedure tvar_keepdims_rank12
        module procedure tvar_keepdims_rank13
        module procedure tvar_keepdims_rank14
        module procedure tvar_keepdims_rank15
    end interface tvar_keepdims

    interface tstd
        module procedure tstd_rank9
        module procedure tstd_rank10
        module procedure tstd_rank11
        module procedure tstd_rank12
        module procedure tstd_rank13
        module procedure tstd_rank14
        module procedure tstd_rank15
    end interface tstd

    interface tstd_all
        module procedure tstd_all_rank9
        module procedure tstd_all_rank10
        module procedure tstd_all_rank11
        module procedure tstd_all_rank12
        module procedure tstd_all_rank13
        module procedure tstd_all_rank14
        module procedure tstd_all_rank15
    end interface tstd_all

    interface tstd_keepdims
        module procedure tstd_keepdims_rank9
        module procedure tstd_keepdims_rank10
        module procedure tstd_keepdims_rank11
        module procedure tstd_keepdims_rank12
        module procedure tstd_keepdims_rank13
        module procedure tstd_keepdims_rank14
        module procedure tstd_keepdims_rank15
    end interface tstd_keepdims

    interface tsem
        module procedure tsem_rank9
        module procedure tsem_rank10
        module procedure tsem_rank11
        module procedure tsem_rank12
        module procedure tsem_rank13
        module procedure tsem_rank14
        module procedure tsem_rank15
    end interface tsem

    interface tsem_all
        module procedure tsem_all_rank9
        module procedure tsem_all_rank10
        module procedure tsem_all_rank11
        module procedure tsem_all_rank12
        module procedure tsem_all_rank13
        module procedure tsem_all_rank14
        module procedure tsem_all_rank15
    end interface tsem_all

    interface tsem_keepdims
        module procedure tsem_keepdims_rank9
        module procedure tsem_keepdims_rank10
        module procedure tsem_keepdims_rank11
        module procedure tsem_keepdims_rank12
        module procedure tsem_keepdims_rank13
        module procedure tsem_keepdims_rank14
        module procedure tsem_keepdims_rank15
    end interface tsem_keepdims

    interface tmin
        module procedure tmin_rank9
        module procedure tmin_rank10
        module procedure tmin_rank11
        module procedure tmin_rank12
        module procedure tmin_rank13
        module procedure tmin_rank14
        module procedure tmin_rank15
    end interface tmin

    interface tmin_all
        module procedure tmin_all_rank9
        module procedure tmin_all_rank10
        module procedure tmin_all_rank11
        module procedure tmin_all_rank12
        module procedure tmin_all_rank13
        module procedure tmin_all_rank14
        module procedure tmin_all_rank15
    end interface tmin_all

    interface tmin_keepdims
        module procedure tmin_keepdims_rank9
        module procedure tmin_keepdims_rank10
        module procedure tmin_keepdims_rank11
        module procedure tmin_keepdims_rank12
        module procedure tmin_keepdims_rank13
        module procedure tmin_keepdims_rank14
        module procedure tmin_keepdims_rank15
    end interface tmin_keepdims

    interface tmax
        module procedure tmax_rank9
        module procedure tmax_rank10
        module procedure tmax_rank11
        module procedure tmax_rank12
        module procedure tmax_rank13
        module procedure tmax_rank14
        module procedure tmax_rank15
    end interface tmax

    interface tmax_all
        module procedure tmax_all_rank9
        module procedure tmax_all_rank10
        module procedure tmax_all_rank11
        module procedure tmax_all_rank12
        module procedure tmax_all_rank13
        module procedure tmax_all_rank14
        module procedure tmax_all_rank15
    end interface tmax_all

    interface tmax_keepdims
        module procedure tmax_keepdims_rank9
        module procedure tmax_keepdims_rank10
        module procedure tmax_keepdims_rank11
        module procedure tmax_keepdims_rank12
        module procedure tmax_keepdims_rank13
        module procedure tmax_keepdims_rank14
        module procedure tmax_keepdims_rank15
    end interface tmax_keepdims

    public :: tmean
    public :: tmean_all
    public :: tmean_keepdims
    public :: tvar
    public :: tvar_all
    public :: tvar_keepdims
    public :: tstd
    public :: tstd_all
    public :: tstd_keepdims
    public :: tsem
    public :: tsem_all
    public :: tsem_keepdims
    public :: tmin
    public :: tmin_all
    public :: tmin_keepdims
    public :: tmax
    public :: tmax_all
    public :: tmax_keepdims

contains

    pure function tmean_all_rank9(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank9

    pure function tmean_axis_rank9(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank9

    pure function tmean_keepdims_rank9(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        real(dp) :: scalar_value

        if (.not. present(axis)) then
            full_shape = 1
            scalar_value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
                upper_limit, inclusive, nan_policy=nan_policy)
            values = reshape([scalar_value], full_shape)
            return
        end if
        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank9

    pure function tvar_all_rank9(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank9

    pure function tvar_rank9(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank9

    pure function tvar_keepdims_rank9(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank9

    pure function tstd_all_rank9(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank9

    pure function tstd_rank9(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank9

    pure function tstd_keepdims_rank9(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank9

    pure function tsem_all_rank9(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank9

    pure function tsem_rank9(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank9

    pure function tsem_keepdims_rank9(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank9

    pure function tmin_all_rank9(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank9

    pure function tmin_rank9(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank9

    pure function tmin_keepdims_rank9(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank9

    pure function tmax_all_rank9(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank9

    pure function tmax_rank9(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank9

    pure function tmax_keepdims_rank9(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank9

    pure function tmean_all_rank10(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank10

    pure function tmean_axis_rank10(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank10

    pure function tmean_keepdims_rank10(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        real(dp) :: scalar_value

        if (.not. present(axis)) then
            full_shape = 1
            scalar_value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
                upper_limit, inclusive, nan_policy=nan_policy)
            values = reshape([scalar_value], full_shape)
            return
        end if
        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank10

    pure function tvar_all_rank10(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank10

    pure function tvar_rank10(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank10

    pure function tvar_keepdims_rank10(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank10

    pure function tstd_all_rank10(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank10

    pure function tstd_rank10(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank10

    pure function tstd_keepdims_rank10(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank10

    pure function tsem_all_rank10(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank10

    pure function tsem_rank10(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank10

    pure function tsem_keepdims_rank10(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank10

    pure function tmin_all_rank10(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank10

    pure function tmin_rank10(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank10

    pure function tmin_keepdims_rank10(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank10

    pure function tmax_all_rank10(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank10

    pure function tmax_rank10(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank10

    pure function tmax_keepdims_rank10(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank10

    pure function tmean_all_rank11(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank11

    pure function tmean_axis_rank11(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank11

    pure function tmean_keepdims_rank11(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        real(dp) :: scalar_value

        if (.not. present(axis)) then
            full_shape = 1
            scalar_value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
                upper_limit, inclusive, nan_policy=nan_policy)
            values = reshape([scalar_value], full_shape)
            return
        end if
        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank11

    pure function tvar_all_rank11(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank11

    pure function tvar_rank11(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank11

    pure function tvar_keepdims_rank11(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank11

    pure function tstd_all_rank11(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank11

    pure function tstd_rank11(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank11

    pure function tstd_keepdims_rank11(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank11

    pure function tsem_all_rank11(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank11

    pure function tsem_rank11(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank11

    pure function tsem_keepdims_rank11(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank11

    pure function tmin_all_rank11(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank11

    pure function tmin_rank11(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank11

    pure function tmin_keepdims_rank11(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank11

    pure function tmax_all_rank11(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank11

    pure function tmax_rank11(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank11

    pure function tmax_keepdims_rank11(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank11

    pure function tmean_all_rank12(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank12

    pure function tmean_axis_rank12(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank12

    pure function tmean_keepdims_rank12(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        real(dp) :: scalar_value

        if (.not. present(axis)) then
            full_shape = 1
            scalar_value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
                upper_limit, inclusive, nan_policy=nan_policy)
            values = reshape([scalar_value], full_shape)
            return
        end if
        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank12

    pure function tvar_all_rank12(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank12

    pure function tvar_rank12(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank12

    pure function tvar_keepdims_rank12(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank12

    pure function tstd_all_rank12(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank12

    pure function tstd_rank12(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank12

    pure function tstd_keepdims_rank12(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank12

    pure function tsem_all_rank12(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank12

    pure function tsem_rank12(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank12

    pure function tsem_keepdims_rank12(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank12

    pure function tmin_all_rank12(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank12

    pure function tmin_rank12(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank12

    pure function tmin_keepdims_rank12(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank12

    pure function tmax_all_rank12(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank12

    pure function tmax_rank12(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank12

    pure function tmax_keepdims_rank12(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank12

    pure function tmean_all_rank13(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank13

    pure function tmean_axis_rank13(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank13

    pure function tmean_keepdims_rank13(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        real(dp) :: scalar_value

        if (.not. present(axis)) then
            full_shape = 1
            scalar_value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
                upper_limit, inclusive, nan_policy=nan_policy)
            values = reshape([scalar_value], full_shape)
            return
        end if
        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank13

    pure function tvar_all_rank13(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank13

    pure function tvar_rank13(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank13

    pure function tvar_keepdims_rank13(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank13

    pure function tstd_all_rank13(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank13

    pure function tstd_rank13(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank13

    pure function tstd_keepdims_rank13(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank13

    pure function tsem_all_rank13(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank13

    pure function tsem_rank13(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank13

    pure function tsem_keepdims_rank13(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank13

    pure function tmin_all_rank13(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank13

    pure function tmin_rank13(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank13

    pure function tmin_keepdims_rank13(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank13

    pure function tmax_all_rank13(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank13

    pure function tmax_rank13(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank13

    pure function tmax_keepdims_rank13(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank13

    pure function tmean_all_rank14(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank14

    pure function tmean_axis_rank14(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank14

    pure function tmean_keepdims_rank14(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        real(dp) :: scalar_value

        if (.not. present(axis)) then
            full_shape = 1
            scalar_value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
                upper_limit, inclusive, nan_policy=nan_policy)
            values = reshape([scalar_value], full_shape)
            return
        end if
        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank14

    pure function tvar_all_rank14(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank14

    pure function tvar_rank14(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank14

    pure function tvar_keepdims_rank14(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank14

    pure function tstd_all_rank14(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank14

    pure function tstd_rank14(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank14

    pure function tstd_keepdims_rank14(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank14

    pure function tsem_all_rank14(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank14

    pure function tsem_rank14(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank14

    pure function tsem_keepdims_rank14(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank14

    pure function tmin_all_rank14(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank14

    pure function tmin_rank14(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank14

    pure function tmin_keepdims_rank14(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank14

    pure function tmax_all_rank14(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank14

    pure function tmax_rank14(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank14

    pure function tmax_keepdims_rank14(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank14

    pure function tmean_all_rank15(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank15

    pure function tmean_axis_rank15(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank15

    pure function tmean_keepdims_rank15(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        real(dp) :: scalar_value

        if (.not. present(axis)) then
            full_shape = 1
            scalar_value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
                upper_limit, inclusive, nan_policy=nan_policy)
            values = reshape([scalar_value], full_shape)
            return
        end if
        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank15

    pure function tvar_all_rank15(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank15

    pure function tvar_rank15(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank15

    pure function tvar_keepdims_rank15(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank15

    pure function tstd_all_rank15(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank15

    pure function tstd_rank15(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank15

    pure function tstd_keepdims_rank15(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank15

    pure function tsem_all_rank15(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank15

    pure function tsem_rank15(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank15

    pure function tsem_keepdims_rank15(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank15

    pure function tmin_all_rank15(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank15

    pure function tmin_rank15(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank15

    pure function tmin_keepdims_rank15(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_min, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank15

    pure function tmax_all_rank15(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank15

    pure function tmax_rank15(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank15

    pure function tmax_keepdims_rank15(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_max, &
            reduced, reduced_shape, ok, &
            upper_limit=upper_limit, inclusive=limit_pair, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank15

end module scifort_descriptive_trimmed_9_15
