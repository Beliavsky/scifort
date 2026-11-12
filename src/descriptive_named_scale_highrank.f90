! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_named_scale_highrank
    use scifort_descriptive_highrank_stats_5_9, only : iqr_5_9 => iqr, &
        mad_5_9 => median_abs_deviation, iqr_keepdims_5_9 => iqr_keepdims, &
        mad_keepdims_5_9 => median_abs_deviation_keepdims
    use scifort_descriptive_highrank_stats_10_15, only : iqr_10_15 => iqr, &
        mad_10_15 => median_abs_deviation, iqr_keepdims_10_15 => iqr_keepdims, &
        mad_keepdims_10_15 => median_abs_deviation_keepdims
    use scifort_descriptive_named_scale_engine, only : named_iqr_scale, named_mad_scale
    use scifort_kinds, only : dp
    implicit none
    private

    interface iqr
        module procedure iqr_named_rank5_all
        module procedure iqr_named_rank5_axis
        module procedure iqr_named_rank6_all
        module procedure iqr_named_rank6_axis
        module procedure iqr_named_rank7_all
        module procedure iqr_named_rank7_axis
        module procedure iqr_named_rank8_all
        module procedure iqr_named_rank8_axis
        module procedure iqr_named_rank9_all
        module procedure iqr_named_rank9_axis
        module procedure iqr_named_rank10_all
        module procedure iqr_named_rank10_axis
        module procedure iqr_named_rank11_all
        module procedure iqr_named_rank11_axis
        module procedure iqr_named_rank12_all
        module procedure iqr_named_rank12_axis
        module procedure iqr_named_rank13_all
        module procedure iqr_named_rank13_axis
        module procedure iqr_named_rank14_all
        module procedure iqr_named_rank14_axis
        module procedure iqr_named_rank15_all
        module procedure iqr_named_rank15_axis
    end interface iqr

    interface median_abs_deviation
        module procedure mad_named_rank5
        module procedure mad_named_rank6
        module procedure mad_named_rank7
        module procedure mad_named_rank8
        module procedure mad_named_rank9
        module procedure mad_named_rank10
        module procedure mad_named_rank11
        module procedure mad_named_rank12
        module procedure mad_named_rank13
        module procedure mad_named_rank14
        module procedure mad_named_rank15
    end interface median_abs_deviation

    interface iqr_keepdims
        module procedure iqr_named_keepdims_rank5
        module procedure iqr_named_keepdims_rank6
        module procedure iqr_named_keepdims_rank7
        module procedure iqr_named_keepdims_rank8
        module procedure iqr_named_keepdims_rank9
        module procedure iqr_named_keepdims_rank10
        module procedure iqr_named_keepdims_rank11
        module procedure iqr_named_keepdims_rank12
        module procedure iqr_named_keepdims_rank13
        module procedure iqr_named_keepdims_rank14
        module procedure iqr_named_keepdims_rank15
    end interface iqr_keepdims

    interface median_abs_deviation_keepdims
        module procedure mad_named_keepdims_rank5
        module procedure mad_named_keepdims_rank6
        module procedure mad_named_keepdims_rank7
        module procedure mad_named_keepdims_rank8
        module procedure mad_named_keepdims_rank9
        module procedure mad_named_keepdims_rank10
        module procedure mad_named_keepdims_rank11
        module procedure mad_named_keepdims_rank12
        module procedure mad_named_keepdims_rank13
        module procedure mad_named_keepdims_rank14
        module procedure mad_named_keepdims_rank15
    end interface median_abs_deviation_keepdims

    public :: iqr
    public :: median_abs_deviation
    public :: iqr_keepdims
    public :: median_abs_deviation_keepdims

contains

    pure function iqr_named_rank5_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_5_9(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank5_all

    pure function iqr_named_rank5_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 5
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank5_axis

    pure function mad_named_rank5(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank5

    pure function iqr_named_keepdims_rank5(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank5

    pure function mad_named_keepdims_rank5(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank5

    pure function iqr_named_rank6_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_5_9(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank6_all

    pure function iqr_named_rank6_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 6
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank6_axis

    pure function mad_named_rank6(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank6

    pure function iqr_named_keepdims_rank6(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank6

    pure function mad_named_keepdims_rank6(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank6

    pure function iqr_named_rank7_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_5_9(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank7_all

    pure function iqr_named_rank7_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 7
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank7_axis

    pure function mad_named_rank7(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank7

    pure function iqr_named_keepdims_rank7(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank7

    pure function mad_named_keepdims_rank7(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank7

    pure function iqr_named_rank8_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_5_9(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank8_all

    pure function iqr_named_rank8_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 8
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank8_axis

    pure function mad_named_rank8(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank8

    pure function iqr_named_keepdims_rank8(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank8

    pure function mad_named_keepdims_rank8(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank8

    pure function iqr_named_rank9_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_5_9(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank9_all

    pure function iqr_named_rank9_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 9
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank9_axis

    pure function mad_named_rank9(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank9

    pure function iqr_named_keepdims_rank9(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_5_9(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank9

    pure function mad_named_keepdims_rank9(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_5_9(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank9

    pure function iqr_named_rank10_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_10_15(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank10_all

    pure function iqr_named_rank10_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 10
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank10_axis

    pure function mad_named_rank10(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank10

    pure function iqr_named_keepdims_rank10(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank10

    pure function mad_named_keepdims_rank10(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank10

    pure function iqr_named_rank11_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_10_15(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank11_all

    pure function iqr_named_rank11_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 11
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank11_axis

    pure function mad_named_rank11(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank11

    pure function iqr_named_keepdims_rank11(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank11

    pure function mad_named_keepdims_rank11(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank11

    pure function iqr_named_rank12_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_10_15(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank12_all

    pure function iqr_named_rank12_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 12
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank12_axis

    pure function mad_named_rank12(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank12

    pure function iqr_named_keepdims_rank12(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank12

    pure function mad_named_keepdims_rank12(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank12

    pure function iqr_named_rank13_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_10_15(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank13_all

    pure function iqr_named_rank13_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 13
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank13_axis

    pure function mad_named_rank13(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank13

    pure function iqr_named_keepdims_rank13(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank13

    pure function mad_named_keepdims_rank13(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank13

    pure function iqr_named_rank14_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_10_15(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank14_all

    pure function iqr_named_rank14_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 14
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank14_axis

    pure function mad_named_rank14(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank14

    pure function iqr_named_keepdims_rank14(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank14

    pure function mad_named_keepdims_rank14(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank14

    pure function iqr_named_rank15_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        value = iqr_10_15(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank15_all

    pure function iqr_named_rank15_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 15
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_rank15_axis

    pure function mad_named_rank15(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_rank15

    pure function iqr_named_keepdims_rank15(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        values = iqr_keepdims_10_15(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
    end function iqr_named_keepdims_rank15

    pure function mad_named_keepdims_rank15(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        values = mad_keepdims_10_15(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
    end function mad_named_keepdims_rank15

end module scifort_descriptive_named_scale_highrank
