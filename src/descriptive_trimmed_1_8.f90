! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_trimmed_1_8
    use scifort_descriptive_trimmed_engine, only : trimmed_kind_max, trimmed_kind_mean, &
        trimmed_kind_min, trimmed_kind_sem, trimmed_kind_std, trimmed_kind_var, &
        trimmed_reduce_assumed_rank, trimmed_stat_all_assumed_rank, trimmed_stat_vector
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface tmean
        module procedure tmean_all_rank1
        module procedure tmean_all_rank2
        module procedure tmean_all_rank3
        module procedure tmean_all_rank4
        module procedure tmean_all_rank5
        module procedure tmean_all_rank6
        module procedure tmean_all_rank7
        module procedure tmean_all_rank8
        module procedure tmean_axis_rank1
        module procedure tmean_axis_rank2
        module procedure tmean_axis_rank3
        module procedure tmean_axis_rank4
        module procedure tmean_axis_rank5
        module procedure tmean_axis_rank6
        module procedure tmean_axis_rank7
        module procedure tmean_axis_rank8
    end interface tmean

    interface tmean_all
        module procedure tmean_all_rank1
        module procedure tmean_all_rank2
        module procedure tmean_all_rank3
        module procedure tmean_all_rank4
        module procedure tmean_all_rank5
        module procedure tmean_all_rank6
        module procedure tmean_all_rank7
        module procedure tmean_all_rank8
    end interface tmean_all

    interface tmean_keepdims
        module procedure tmean_keepdims_rank1
        module procedure tmean_keepdims_rank2
        module procedure tmean_keepdims_rank3
        module procedure tmean_keepdims_rank4
        module procedure tmean_keepdims_rank5
        module procedure tmean_keepdims_rank6
        module procedure tmean_keepdims_rank7
        module procedure tmean_keepdims_rank8
    end interface tmean_keepdims

    interface tvar
        module procedure tvar_rank1
        module procedure tvar_rank2
        module procedure tvar_rank3
        module procedure tvar_rank4
        module procedure tvar_rank5
        module procedure tvar_rank6
        module procedure tvar_rank7
        module procedure tvar_rank8
    end interface tvar

    interface tvar_all
        module procedure tvar_all_rank1
        module procedure tvar_all_rank2
        module procedure tvar_all_rank3
        module procedure tvar_all_rank4
        module procedure tvar_all_rank5
        module procedure tvar_all_rank6
        module procedure tvar_all_rank7
        module procedure tvar_all_rank8
    end interface tvar_all

    interface tvar_keepdims
        module procedure tvar_keepdims_rank1
        module procedure tvar_keepdims_rank2
        module procedure tvar_keepdims_rank3
        module procedure tvar_keepdims_rank4
        module procedure tvar_keepdims_rank5
        module procedure tvar_keepdims_rank6
        module procedure tvar_keepdims_rank7
        module procedure tvar_keepdims_rank8
    end interface tvar_keepdims

    interface tstd
        module procedure tstd_rank1
        module procedure tstd_rank2
        module procedure tstd_rank3
        module procedure tstd_rank4
        module procedure tstd_rank5
        module procedure tstd_rank6
        module procedure tstd_rank7
        module procedure tstd_rank8
    end interface tstd

    interface tstd_all
        module procedure tstd_all_rank1
        module procedure tstd_all_rank2
        module procedure tstd_all_rank3
        module procedure tstd_all_rank4
        module procedure tstd_all_rank5
        module procedure tstd_all_rank6
        module procedure tstd_all_rank7
        module procedure tstd_all_rank8
    end interface tstd_all

    interface tstd_keepdims
        module procedure tstd_keepdims_rank1
        module procedure tstd_keepdims_rank2
        module procedure tstd_keepdims_rank3
        module procedure tstd_keepdims_rank4
        module procedure tstd_keepdims_rank5
        module procedure tstd_keepdims_rank6
        module procedure tstd_keepdims_rank7
        module procedure tstd_keepdims_rank8
    end interface tstd_keepdims

    interface tsem
        module procedure tsem_rank1
        module procedure tsem_rank2
        module procedure tsem_rank3
        module procedure tsem_rank4
        module procedure tsem_rank5
        module procedure tsem_rank6
        module procedure tsem_rank7
        module procedure tsem_rank8
    end interface tsem

    interface tsem_all
        module procedure tsem_all_rank1
        module procedure tsem_all_rank2
        module procedure tsem_all_rank3
        module procedure tsem_all_rank4
        module procedure tsem_all_rank5
        module procedure tsem_all_rank6
        module procedure tsem_all_rank7
        module procedure tsem_all_rank8
    end interface tsem_all

    interface tsem_keepdims
        module procedure tsem_keepdims_rank1
        module procedure tsem_keepdims_rank2
        module procedure tsem_keepdims_rank3
        module procedure tsem_keepdims_rank4
        module procedure tsem_keepdims_rank5
        module procedure tsem_keepdims_rank6
        module procedure tsem_keepdims_rank7
        module procedure tsem_keepdims_rank8
    end interface tsem_keepdims

    interface tmin
        module procedure tmin_rank1
        module procedure tmin_rank2
        module procedure tmin_rank3
        module procedure tmin_rank4
        module procedure tmin_rank5
        module procedure tmin_rank6
        module procedure tmin_rank7
        module procedure tmin_rank8
    end interface tmin

    interface tmin_all
        module procedure tmin_all_rank1
        module procedure tmin_all_rank2
        module procedure tmin_all_rank3
        module procedure tmin_all_rank4
        module procedure tmin_all_rank5
        module procedure tmin_all_rank6
        module procedure tmin_all_rank7
        module procedure tmin_all_rank8
    end interface tmin_all

    interface tmin_keepdims
        module procedure tmin_keepdims_rank1
        module procedure tmin_keepdims_rank2
        module procedure tmin_keepdims_rank3
        module procedure tmin_keepdims_rank4
        module procedure tmin_keepdims_rank5
        module procedure tmin_keepdims_rank6
        module procedure tmin_keepdims_rank7
        module procedure tmin_keepdims_rank8
    end interface tmin_keepdims

    interface tmax
        module procedure tmax_rank1
        module procedure tmax_rank2
        module procedure tmax_rank3
        module procedure tmax_rank4
        module procedure tmax_rank5
        module procedure tmax_rank6
        module procedure tmax_rank7
        module procedure tmax_rank8
    end interface tmax

    interface tmax_all
        module procedure tmax_all_rank1
        module procedure tmax_all_rank2
        module procedure tmax_all_rank3
        module procedure tmax_all_rank4
        module procedure tmax_all_rank5
        module procedure tmax_all_rank6
        module procedure tmax_all_rank7
        module procedure tmax_all_rank8
    end interface tmax_all

    interface tmax_keepdims
        module procedure tmax_keepdims_rank1
        module procedure tmax_keepdims_rank2
        module procedure tmax_keepdims_rank3
        module procedure tmax_keepdims_rank4
        module procedure tmax_keepdims_rank5
        module procedure tmax_keepdims_rank6
        module procedure tmax_keepdims_rank7
        module procedure tmax_keepdims_rank8
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

    pure function tmean_all_rank1(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank1

    pure function tmean_axis_rank1(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: axis_value

        axis_value = axis
        if (axis_value /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = trimmed_stat_vector(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_axis_rank1

    pure function tmean_keepdims_rank1(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        real(dp) :: scalar_value

        if (.not. present(axis)) then
            full_shape = 1
            scalar_value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
                upper_limit, inclusive, nan_policy=nan_policy)
            values = reshape([scalar_value], full_shape)
            return
        end if
        axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        full_shape = 1
        scalar_value = trimmed_stat_vector(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
        values = reshape([scalar_value], full_shape)
    end function tmean_keepdims_rank1

    pure function tvar_all_rank1(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank1

    pure function tvar_rank1(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = trimmed_stat_vector(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_rank1

    pure function tvar_keepdims_rank1(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        real(dp) :: scalar_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        full_shape = 1
        scalar_value = trimmed_stat_vector(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
        values = reshape([scalar_value], full_shape)
    end function tvar_keepdims_rank1

    pure function tstd_all_rank1(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank1

    pure function tstd_rank1(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = trimmed_stat_vector(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_rank1

    pure function tstd_keepdims_rank1(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        real(dp) :: scalar_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        full_shape = 1
        scalar_value = trimmed_stat_vector(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
        values = reshape([scalar_value], full_shape)
    end function tstd_keepdims_rank1

    pure function tsem_all_rank1(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank1

    pure function tsem_rank1(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = trimmed_stat_vector(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_rank1

    pure function tsem_keepdims_rank1(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        real(dp) :: scalar_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        full_shape = 1
        scalar_value = trimmed_stat_vector(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
        values = reshape([scalar_value], full_shape)
    end function tsem_keepdims_rank1

    pure function tmin_all_rank1(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank1

    pure function tmin_rank1(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: axis_value
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        if (axis_value /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = trimmed_stat_vector(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_rank1

    pure function tmin_keepdims_rank1(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        logical :: limit_pair(2)
        real(dp) :: scalar_value

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        full_shape = 1
        scalar_value = trimmed_stat_vector(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
        values = reshape([scalar_value], full_shape)
    end function tmin_keepdims_rank1

    pure function tmax_all_rank1(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank1

    pure function tmax_rank1(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: axis_value
        logical :: limit_pair(2)

        axis_value = 1
        if (present(axis)) axis_value = axis
        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        if (axis_value /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = trimmed_stat_vector(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_rank1

    pure function tmax_keepdims_rank1(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        logical :: limit_pair(2)
        real(dp) :: scalar_value

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        full_shape = 1
        scalar_value = trimmed_stat_vector(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
        values = reshape([scalar_value], full_shape)
    end function tmax_keepdims_rank1

    pure function tmean_all_rank2(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank2

    pure function tmean_axis_rank2(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:)

        integer :: axis_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank2

    pure function tmean_keepdims_rank2(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
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
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank2

    pure function tvar_all_rank2(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank2

    pure function tvar_rank2(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:)

        integer :: axis_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank2

    pure function tvar_keepdims_rank2(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank2

    pure function tstd_all_rank2(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank2

    pure function tstd_rank2(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:)

        integer :: axis_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank2

    pure function tstd_keepdims_rank2(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank2

    pure function tsem_all_rank2(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank2

    pure function tsem_rank2(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:)

        integer :: axis_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank2

    pure function tsem_keepdims_rank2(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank2

    pure function tmin_all_rank2(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank2

    pure function tmin_rank2(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:)

        integer :: axis_value
        integer :: reduced_shape(1)
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
            allocate(value(0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank2

    pure function tmin_keepdims_rank2(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
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
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank2

    pure function tmax_all_rank2(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank2

    pure function tmax_rank2(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:)

        integer :: axis_value
        integer :: reduced_shape(1)
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
            allocate(value(0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank2

    pure function tmax_keepdims_rank2(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
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
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank2

    pure function tmean_all_rank3(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank3

    pure function tmean_axis_rank3(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :)

        integer :: axis_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank3

    pure function tmean_keepdims_rank3(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
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
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank3

    pure function tvar_all_rank3(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank3

    pure function tvar_rank3(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :)

        integer :: axis_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank3

    pure function tvar_keepdims_rank3(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank3

    pure function tstd_all_rank3(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank3

    pure function tstd_rank3(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :)

        integer :: axis_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank3

    pure function tstd_keepdims_rank3(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank3

    pure function tsem_all_rank3(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank3

    pure function tsem_rank3(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :)

        integer :: axis_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank3

    pure function tsem_keepdims_rank3(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank3

    pure function tmin_all_rank3(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank3

    pure function tmin_rank3(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :)

        integer :: axis_value
        integer :: reduced_shape(2)
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
            allocate(value(0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank3

    pure function tmin_keepdims_rank3(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
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
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank3

    pure function tmax_all_rank3(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank3

    pure function tmax_rank3(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :)

        integer :: axis_value
        integer :: reduced_shape(2)
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
            allocate(value(0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank3

    pure function tmax_keepdims_rank3(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
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
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank3

    pure function tmean_all_rank4(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank4

    pure function tmean_axis_rank4(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :)

        integer :: axis_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank4

    pure function tmean_keepdims_rank4(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
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
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank4

    pure function tvar_all_rank4(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank4

    pure function tvar_rank4(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :)

        integer :: axis_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank4

    pure function tvar_keepdims_rank4(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank4

    pure function tstd_all_rank4(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank4

    pure function tstd_rank4(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :)

        integer :: axis_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank4

    pure function tstd_keepdims_rank4(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank4

    pure function tsem_all_rank4(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank4

    pure function tsem_rank4(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :)

        integer :: axis_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank4

    pure function tsem_keepdims_rank4(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank4

    pure function tmin_all_rank4(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank4

    pure function tmin_rank4(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :)

        integer :: axis_value
        integer :: reduced_shape(3)
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
            allocate(value(0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank4

    pure function tmin_keepdims_rank4(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
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
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank4

    pure function tmax_all_rank4(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank4

    pure function tmax_rank4(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :)

        integer :: axis_value
        integer :: reduced_shape(3)
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
            allocate(value(0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank4

    pure function tmax_keepdims_rank4(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
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
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank4

    pure function tmean_all_rank5(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank5

    pure function tmean_axis_rank5(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank5

    pure function tmean_keepdims_rank5(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
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
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank5

    pure function tvar_all_rank5(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank5

    pure function tvar_rank5(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank5

    pure function tvar_keepdims_rank5(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank5

    pure function tstd_all_rank5(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank5

    pure function tstd_rank5(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank5

    pure function tstd_keepdims_rank5(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank5

    pure function tsem_all_rank5(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank5

    pure function tsem_rank5(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank5

    pure function tsem_keepdims_rank5(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank5

    pure function tmin_all_rank5(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank5

    pure function tmin_rank5(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
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
            allocate(value(0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank5

    pure function tmin_keepdims_rank5(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
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
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank5

    pure function tmax_all_rank5(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank5

    pure function tmax_rank5(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
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
            allocate(value(0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank5

    pure function tmax_keepdims_rank5(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
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
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank5

    pure function tmean_all_rank6(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank6

    pure function tmean_axis_rank6(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank6

    pure function tmean_keepdims_rank6(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
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
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank6

    pure function tvar_all_rank6(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank6

    pure function tvar_rank6(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank6

    pure function tvar_keepdims_rank6(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank6

    pure function tstd_all_rank6(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank6

    pure function tstd_rank6(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank6

    pure function tstd_keepdims_rank6(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank6

    pure function tsem_all_rank6(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank6

    pure function tsem_rank6(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank6

    pure function tsem_keepdims_rank6(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank6

    pure function tmin_all_rank6(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank6

    pure function tmin_rank6(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
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
            allocate(value(0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank6

    pure function tmin_keepdims_rank6(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
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
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank6

    pure function tmax_all_rank6(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank6

    pure function tmax_rank6(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
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
            allocate(value(0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank6

    pure function tmax_keepdims_rank6(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
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
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank6

    pure function tmean_all_rank7(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank7

    pure function tmean_axis_rank7(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank7

    pure function tmean_keepdims_rank7(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
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
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank7

    pure function tvar_all_rank7(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank7

    pure function tvar_rank7(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank7

    pure function tvar_keepdims_rank7(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank7

    pure function tstd_all_rank7(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank7

    pure function tstd_rank7(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank7

    pure function tstd_keepdims_rank7(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank7

    pure function tsem_all_rank7(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank7

    pure function tsem_rank7(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank7

    pure function tsem_keepdims_rank7(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank7

    pure function tmin_all_rank7(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank7

    pure function tmin_rank7(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
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
            allocate(value(0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank7

    pure function tmin_keepdims_rank7(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
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
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank7

    pure function tmax_all_rank7(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank7

    pure function tmax_rank7(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
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
            allocate(value(0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank7

    pure function tmax_keepdims_rank7(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
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
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank7

    pure function tmean_all_rank8(x, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_mean, lower_limit, &
            upper_limit, inclusive, nan_policy=nan_policy)
    end function tmean_all_rank8

    pure function tmean_axis_rank8(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_mean, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmean_axis_rank8

    pure function tmean_keepdims_rank8(x, axis, lower_limit, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
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
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmean_keepdims_rank8

    pure function tvar_all_rank8(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_var, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tvar_all_rank8

    pure function tvar_rank8(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tvar_rank8

    pure function tvar_keepdims_rank8(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_var, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tvar_keepdims_rank8

    pure function tstd_all_rank8(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_std, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tstd_all_rank8

    pure function tstd_rank8(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tstd_rank8

    pure function tstd_keepdims_rank8(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_std, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tstd_keepdims_rank8

    pure function tsem_all_rank8(x, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_sem, lower_limit, &
            upper_limit, inclusive, ddof=ddof, &
            nan_policy=nan_policy)
    end function tsem_all_rank8

    pure function tsem_rank8(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tsem_rank8

    pure function tsem_keepdims_rank8(x, axis, lower_limit, upper_limit, inclusive, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive(2) !! include exact lower/upper bounds
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trimmed_reduce_assumed_rank(x, axis_value, trimmed_kind_sem, &
            reduced, reduced_shape, ok, &
            lower_limit=lower_limit, upper_limit=upper_limit, inclusive=inclusive, &
            ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tsem_keepdims_rank8

    pure function tmin_all_rank8(x, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(1) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_min, lower_limit=lower_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmin_all_rank8

    pure function tmin_rank8(x, axis, lower_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
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
            allocate(value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmin_rank8

    pure function tmin_keepdims_rank8(x, axis, lower_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: lower_limit !! optional lower value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
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
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmin_keepdims_rank8

    pure function tmax_all_rank8(x, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        logical :: limit_pair(2)

        limit_pair = .true.
        if (present(inclusive)) limit_pair(2) = inclusive
        value = trimmed_stat_all_assumed_rank(x, trimmed_kind_max, upper_limit=upper_limit, &
            inclusive=limit_pair, nan_policy=nan_policy)
    end function tmax_all_rank8

    pure function tmax_rank8(x, axis, upper_limit, inclusive, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
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
            allocate(value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function tmax_rank8

    pure function tmax_keepdims_rank8(x, axis, upper_limit, inclusive, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        real(dp), intent(in), optional :: upper_limit !! optional upper value bound
        logical, intent(in), optional :: inclusive !! include exact limit value
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
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
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function tmax_keepdims_rank8

end module scifort_descriptive_trimmed_1_8
