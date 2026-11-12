! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_named_scale_lowrank
    use scifort_descriptive_lowrank, only : iqr_numeric => iqr, &
        mad_numeric => median_abs_deviation, iqr_keepdims_numeric => iqr_keepdims, &
        mad_keepdims_numeric => median_abs_deviation_keepdims
    use scifort_descriptive_named_scale_engine, only : named_iqr_scale, named_mad_scale
    use scifort_kinds, only : dp
    implicit none
    private

    interface iqr
        module procedure iqr_named_rank1_all
        module procedure iqr_named_rank2_all
        module procedure iqr_named_rank2_axis
        module procedure iqr_named_rank3_all
        module procedure iqr_named_rank3_axis
        module procedure iqr_named_rank4_all
        module procedure iqr_named_rank4_axis
    end interface iqr

    interface median_abs_deviation
        module procedure mad_named_rank1
        module procedure mad_named_rank2
        module procedure mad_named_rank3
        module procedure mad_named_rank4
    end interface median_abs_deviation

    interface iqr_keepdims
        module procedure iqr_named_keepdims_rank1
        module procedure iqr_named_keepdims_rank2
        module procedure iqr_named_keepdims_rank3
        module procedure iqr_named_keepdims_rank4
    end interface iqr_keepdims

    interface median_abs_deviation_keepdims
        module procedure mad_named_keepdims_rank1
        module procedure mad_named_keepdims_rank2
        module procedure mad_named_keepdims_rank3
        module procedure mad_named_keepdims_rank4
    end interface median_abs_deviation_keepdims

    public :: iqr
    public :: median_abs_deviation
    public :: iqr_keepdims
    public :: median_abs_deviation_keepdims

contains

    pure function iqr_named_rank1_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            value = iqr_numeric(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            value = iqr_numeric(x, rng=rng, scale=scale_value)
        end if
    end function iqr_named_rank1_all

    pure function mad_named_rank1(x, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        if (present(nan_policy)) then
            value = mad_numeric(x, scale=scale_value, nan_policy=nan_policy)
        else
            value = mad_numeric(x, scale=scale_value)
        end if
    end function mad_named_rank1

    pure function iqr_named_keepdims_rank1(x, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            values = iqr_keepdims_numeric(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            values = iqr_keepdims_numeric(x, rng=rng, scale=scale_value)
        end if
    end function iqr_named_keepdims_rank1

    pure function mad_named_keepdims_rank1(x, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        if (present(nan_policy)) then
            values = mad_keepdims_numeric(x, scale=scale_value, nan_policy=nan_policy)
        else
            values = mad_keepdims_numeric(x, scale=scale_value)
        end if
    end function mad_named_keepdims_rank1

    pure function iqr_named_rank2_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            value = iqr_numeric(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            value = iqr_numeric(x, rng=rng, scale=scale_value)
        end if
    end function iqr_named_rank2_all

    pure function iqr_named_rank2_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 2
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            values = iqr_numeric(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            values = iqr_numeric(x, axis=axis, rng=rng, scale=scale_value)
        end if
    end function iqr_named_rank2_axis

    pure function mad_named_rank2(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 2; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        if (present(nan_policy)) then
            values = mad_numeric(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
        else
            values = mad_numeric(x, axis=axis, scale=scale_value)
        end if
    end function mad_named_rank2

    pure function iqr_named_keepdims_rank2(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 2; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            values = iqr_keepdims_numeric(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            values = iqr_keepdims_numeric(x, axis=axis, rng=rng, scale=scale_value)
        end if
    end function iqr_named_keepdims_rank2

    pure function mad_named_keepdims_rank2(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 2; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        if (present(nan_policy)) then
            values = mad_keepdims_numeric(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
        else
            values = mad_keepdims_numeric(x, axis=axis, scale=scale_value)
        end if
    end function mad_named_keepdims_rank2

    pure function iqr_named_rank3_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            value = iqr_numeric(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            value = iqr_numeric(x, rng=rng, scale=scale_value)
        end if
    end function iqr_named_rank3_all

    pure function iqr_named_rank3_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 3
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            values = iqr_numeric(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            values = iqr_numeric(x, axis=axis, rng=rng, scale=scale_value)
        end if
    end function iqr_named_rank3_axis

    pure function mad_named_rank3(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        if (present(nan_policy)) then
            values = mad_numeric(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
        else
            values = mad_numeric(x, axis=axis, scale=scale_value)
        end if
    end function mad_named_rank3

    pure function iqr_named_keepdims_rank3(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            values = iqr_keepdims_numeric(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            values = iqr_keepdims_numeric(x, axis=axis, rng=rng, scale=scale_value)
        end if
    end function iqr_named_keepdims_rank3

    pure function mad_named_keepdims_rank3(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        if (present(nan_policy)) then
            values = mad_keepdims_numeric(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
        else
            values = mad_keepdims_numeric(x, axis=axis, scale=scale_value)
        end if
    end function mad_named_keepdims_rank3

    pure function iqr_named_rank4_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            value = iqr_numeric(x, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            value = iqr_numeric(x, rng=rng, scale=scale_value)
        end if
    end function iqr_named_rank4_all

    pure function iqr_named_rank4_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in) :: axis !! Fortran axis, 1 to 4
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            values = iqr_numeric(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            values = iqr_numeric(x, axis=axis, rng=rng, scale=scale_value)
        end if
    end function iqr_named_rank4_axis

    pure function mad_named_rank4(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        if (present(nan_policy)) then
            values = mad_numeric(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
        else
            values = mad_numeric(x, axis=axis, scale=scale_value)
        end if
    end function mad_named_rank4

    pure function iqr_named_keepdims_rank4(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp) :: scale_value

        scale_value = named_iqr_scale(scale)
        if (present(nan_policy)) then
            values = iqr_keepdims_numeric(x, axis=axis, rng=rng, scale=scale_value, nan_policy=nan_policy)
        else
            values = iqr_keepdims_numeric(x, axis=axis, rng=rng, scale=scale_value)
        end if
    end function iqr_named_keepdims_rank4

    pure function mad_named_keepdims_rank4(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        character(len=*), intent(in) :: scale !! named divisor; currently 'normal'
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)
        real(dp) :: scale_value

        scale_value = named_mad_scale(scale)
        if (present(nan_policy)) then
            values = mad_keepdims_numeric(x, axis=axis, scale=scale_value, nan_policy=nan_policy)
        else
            values = mad_keepdims_numeric(x, axis=axis, scale=scale_value)
        end if
    end function mad_named_keepdims_rank4

end module scifort_descriptive_named_scale_lowrank
