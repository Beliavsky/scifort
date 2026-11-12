! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_highrank_stats_5_9
    use scifort_descriptive_highrank_stats_engine
    use scifort_kinds, only : dp
    implicit none
    private

    interface sem
        module procedure sem_rank5_assumed
        module procedure sem_rank6_assumed
        module procedure sem_rank7_assumed
        module procedure sem_rank8_assumed
        module procedure sem_rank9_assumed
    end interface sem

    interface variation
        module procedure variation_rank5_assumed
        module procedure variation_rank6_assumed
        module procedure variation_rank7_assumed
        module procedure variation_rank8_assumed
        module procedure variation_rank9_assumed
    end interface variation

    interface median_abs_deviation
        module procedure median_abs_deviation_rank5_assumed
        module procedure median_abs_deviation_rank6_assumed
        module procedure median_abs_deviation_rank7_assumed
        module procedure median_abs_deviation_rank8_assumed
        module procedure median_abs_deviation_rank9_assumed
    end interface median_abs_deviation

    interface moment
        module procedure moment_rank5_assumed
        module procedure moment_rank6_assumed
        module procedure moment_rank7_assumed
        module procedure moment_rank8_assumed
        module procedure moment_rank9_assumed
    end interface moment

    interface skew
        module procedure skew_rank5_assumed
        module procedure skew_rank6_assumed
        module procedure skew_rank7_assumed
        module procedure skew_rank8_assumed
        module procedure skew_rank9_assumed
    end interface skew

    interface kurtosis
        module procedure kurtosis_rank5_assumed
        module procedure kurtosis_rank6_assumed
        module procedure kurtosis_rank7_assumed
        module procedure kurtosis_rank8_assumed
        module procedure kurtosis_rank9_assumed
    end interface kurtosis

    interface iqr
        module procedure iqr_rank5_assumed_all
        module procedure iqr_rank5_assumed_axis
        module procedure iqr_rank6_assumed_all
        module procedure iqr_rank6_assumed_axis
        module procedure iqr_rank7_assumed_all
        module procedure iqr_rank7_assumed_axis
        module procedure iqr_rank8_assumed_all
        module procedure iqr_rank8_assumed_axis
        module procedure iqr_rank9_assumed_all
        module procedure iqr_rank9_assumed_axis
    end interface iqr

    interface sem_keepdims
        module procedure sem_keepdims_rank5_assumed
        module procedure sem_keepdims_rank6_assumed
        module procedure sem_keepdims_rank7_assumed
        module procedure sem_keepdims_rank8_assumed
        module procedure sem_keepdims_rank9_assumed
    end interface sem_keepdims

    interface variation_keepdims
        module procedure variation_keepdims_rank5_assumed
        module procedure variation_keepdims_rank6_assumed
        module procedure variation_keepdims_rank7_assumed
        module procedure variation_keepdims_rank8_assumed
        module procedure variation_keepdims_rank9_assumed
    end interface variation_keepdims

    interface median_abs_deviation_keepdims
        module procedure median_abs_deviation_keepdims_rank5_assumed
        module procedure median_abs_deviation_keepdims_rank6_assumed
        module procedure median_abs_deviation_keepdims_rank7_assumed
        module procedure median_abs_deviation_keepdims_rank8_assumed
        module procedure median_abs_deviation_keepdims_rank9_assumed
    end interface median_abs_deviation_keepdims

    interface moment_keepdims
        module procedure moment_keepdims_rank5_assumed
        module procedure moment_keepdims_rank6_assumed
        module procedure moment_keepdims_rank7_assumed
        module procedure moment_keepdims_rank8_assumed
        module procedure moment_keepdims_rank9_assumed
    end interface moment_keepdims

    interface skew_keepdims
        module procedure skew_keepdims_rank5_assumed
        module procedure skew_keepdims_rank6_assumed
        module procedure skew_keepdims_rank7_assumed
        module procedure skew_keepdims_rank8_assumed
        module procedure skew_keepdims_rank9_assumed
    end interface skew_keepdims

    interface kurtosis_keepdims
        module procedure kurtosis_keepdims_rank5_assumed
        module procedure kurtosis_keepdims_rank6_assumed
        module procedure kurtosis_keepdims_rank7_assumed
        module procedure kurtosis_keepdims_rank8_assumed
        module procedure kurtosis_keepdims_rank9_assumed
    end interface kurtosis_keepdims

    interface iqr_keepdims
        module procedure iqr_keepdims_rank5_assumed
        module procedure iqr_keepdims_rank6_assumed
        module procedure iqr_keepdims_rank7_assumed
        module procedure iqr_keepdims_rank8_assumed
        module procedure iqr_keepdims_rank9_assumed
    end interface iqr_keepdims

    public :: sem
    public :: variation
    public :: median_abs_deviation
    public :: moment
    public :: skew
    public :: kurtosis
    public :: iqr
    public :: sem_keepdims
    public :: variation_keepdims
    public :: median_abs_deviation_keepdims
    public :: moment_keepdims
    public :: skew_keepdims
    public :: kurtosis_keepdims
    public :: iqr_keepdims

contains

    pure function sem_rank5_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function sem_rank5_assumed

    pure function variation_rank5_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function variation_rank5_assumed

    pure function median_abs_deviation_rank5_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function median_abs_deviation_rank5_assumed

    pure function moment_rank5_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function moment_rank5_assumed

    pure function skew_rank5_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function skew_rank5_assumed

    pure function kurtosis_rank5_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kurtosis_rank5_assumed

    pure function iqr_rank5_assumed_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = iqr_assumed_rank_all(x, rng=rng, scale=scale, nan_policy=nan_policy)
    end function iqr_rank5_assumed_all

    pure function iqr_rank5_assumed_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 5
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call reduce_assumed_rank_stat(x, axis, stat_kind_iqr, reduced, reduced_shape, &
            ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function iqr_rank5_assumed_axis

    pure function sem_keepdims_rank5_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function sem_keepdims_rank5_assumed

    pure function variation_keepdims_rank5_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function variation_keepdims_rank5_assumed

    pure function median_abs_deviation_keepdims_rank5_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function median_abs_deviation_keepdims_rank5_assumed

    pure function moment_keepdims_rank5_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function moment_keepdims_rank5_assumed

    pure function skew_keepdims_rank5_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function skew_keepdims_rank5_assumed

    pure function kurtosis_keepdims_rank5_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function kurtosis_keepdims_rank5_assumed

    pure function iqr_keepdims_rank5_assumed(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_iqr, reduced, &
            reduced_shape, ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function iqr_keepdims_rank5_assumed

    pure function sem_rank6_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function sem_rank6_assumed

    pure function variation_rank6_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function variation_rank6_assumed

    pure function median_abs_deviation_rank6_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function median_abs_deviation_rank6_assumed

    pure function moment_rank6_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function moment_rank6_assumed

    pure function skew_rank6_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function skew_rank6_assumed

    pure function kurtosis_rank6_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kurtosis_rank6_assumed

    pure function iqr_rank6_assumed_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = iqr_assumed_rank_all(x, rng=rng, scale=scale, nan_policy=nan_policy)
    end function iqr_rank6_assumed_all

    pure function iqr_rank6_assumed_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 6
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call reduce_assumed_rank_stat(x, axis, stat_kind_iqr, reduced, reduced_shape, &
            ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function iqr_rank6_assumed_axis

    pure function sem_keepdims_rank6_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function sem_keepdims_rank6_assumed

    pure function variation_keepdims_rank6_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function variation_keepdims_rank6_assumed

    pure function median_abs_deviation_keepdims_rank6_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function median_abs_deviation_keepdims_rank6_assumed

    pure function moment_keepdims_rank6_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function moment_keepdims_rank6_assumed

    pure function skew_keepdims_rank6_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function skew_keepdims_rank6_assumed

    pure function kurtosis_keepdims_rank6_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function kurtosis_keepdims_rank6_assumed

    pure function iqr_keepdims_rank6_assumed(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_iqr, reduced, &
            reduced_shape, ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function iqr_keepdims_rank6_assumed

    pure function sem_rank7_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function sem_rank7_assumed

    pure function variation_rank7_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function variation_rank7_assumed

    pure function median_abs_deviation_rank7_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function median_abs_deviation_rank7_assumed

    pure function moment_rank7_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function moment_rank7_assumed

    pure function skew_rank7_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function skew_rank7_assumed

    pure function kurtosis_rank7_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kurtosis_rank7_assumed

    pure function iqr_rank7_assumed_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = iqr_assumed_rank_all(x, rng=rng, scale=scale, nan_policy=nan_policy)
    end function iqr_rank7_assumed_all

    pure function iqr_rank7_assumed_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 7
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call reduce_assumed_rank_stat(x, axis, stat_kind_iqr, reduced, reduced_shape, &
            ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function iqr_rank7_assumed_axis

    pure function sem_keepdims_rank7_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function sem_keepdims_rank7_assumed

    pure function variation_keepdims_rank7_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function variation_keepdims_rank7_assumed

    pure function median_abs_deviation_keepdims_rank7_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function median_abs_deviation_keepdims_rank7_assumed

    pure function moment_keepdims_rank7_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function moment_keepdims_rank7_assumed

    pure function skew_keepdims_rank7_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function skew_keepdims_rank7_assumed

    pure function kurtosis_keepdims_rank7_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function kurtosis_keepdims_rank7_assumed

    pure function iqr_keepdims_rank7_assumed(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_iqr, reduced, &
            reduced_shape, ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function iqr_keepdims_rank7_assumed

    pure function sem_rank8_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function sem_rank8_assumed

    pure function variation_rank8_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function variation_rank8_assumed

    pure function median_abs_deviation_rank8_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function median_abs_deviation_rank8_assumed

    pure function moment_rank8_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function moment_rank8_assumed

    pure function skew_rank8_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function skew_rank8_assumed

    pure function kurtosis_rank8_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kurtosis_rank8_assumed

    pure function iqr_rank8_assumed_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = iqr_assumed_rank_all(x, rng=rng, scale=scale, nan_policy=nan_policy)
    end function iqr_rank8_assumed_all

    pure function iqr_rank8_assumed_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 8
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call reduce_assumed_rank_stat(x, axis, stat_kind_iqr, reduced, reduced_shape, &
            ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function iqr_rank8_assumed_axis

    pure function sem_keepdims_rank8_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function sem_keepdims_rank8_assumed

    pure function variation_keepdims_rank8_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function variation_keepdims_rank8_assumed

    pure function median_abs_deviation_keepdims_rank8_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function median_abs_deviation_keepdims_rank8_assumed

    pure function moment_keepdims_rank8_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function moment_keepdims_rank8_assumed

    pure function skew_keepdims_rank8_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function skew_keepdims_rank8_assumed

    pure function kurtosis_keepdims_rank8_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function kurtosis_keepdims_rank8_assumed

    pure function iqr_keepdims_rank8_assumed(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_iqr, reduced, &
            reduced_shape, ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function iqr_keepdims_rank8_assumed

    pure function sem_rank9_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function sem_rank9_assumed

    pure function variation_rank9_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function variation_rank9_assumed

    pure function median_abs_deviation_rank9_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function median_abs_deviation_rank9_assumed

    pure function moment_rank9_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function moment_rank9_assumed

    pure function skew_rank9_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function skew_rank9_assumed

    pure function kurtosis_rank9_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function kurtosis_rank9_assumed

    pure function iqr_rank9_assumed_all(x, rng, scale, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations reduced over all axes
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = iqr_assumed_rank_all(x, rng=rng, scale=scale, nan_policy=nan_policy)
    end function iqr_rank9_assumed_all

    pure function iqr_rank9_assumed_axis(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced along axis
        integer, intent(in) :: axis !! Fortran axis, 1 to 9
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call reduce_assumed_rank_stat(x, axis, stat_kind_iqr, reduced, reduced_shape, &
            ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function iqr_rank9_assumed_axis

    pure function sem_keepdims_rank9_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_sem, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function sem_keepdims_rank9_assumed

    pure function variation_keepdims_rank9_assumed(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_variation, reduced, &
            reduced_shape, ok, ddof=ddof, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function variation_keepdims_rank9_assumed

    pure function median_abs_deviation_keepdims_rank9_assumed(x, axis, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in), optional :: scale !! divisor applied to MAD, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_mad, reduced, &
            reduced_shape, ok, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function median_abs_deviation_keepdims_rank9_assumed

    pure function moment_keepdims_rank9_assumed(x, order, axis, center, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: order !! nonnegative moment order, default 1
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in), optional :: center !! finite center; default slice mean
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_moment, reduced, &
            reduced_shape, ok, order=order, center=center, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function moment_keepdims_rank9_assumed

    pure function skew_keepdims_rank9_assumed(x, axis, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_skew, reduced, &
            reduced_shape, ok, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function skew_keepdims_rank9_assumed

    pure function kurtosis_keepdims_rank9_assumed(x, axis, fisher, bias, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        logical, intent(in), optional :: fisher !! subtract 3, default true
        logical, intent(in), optional :: bias !! use biased sample moments, default true
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_kurtosis, reduced, &
            reduced_shape, ok, fisher=fisher, bias=bias, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function kurtosis_keepdims_rank9_assumed

    pure function iqr_keepdims_rank9_assumed(x, axis, rng, scale, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; reduced with dimension retained
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in), optional :: rng(2) !! percentile endpoints, default [25,75]
        real(dp), intent(in), optional :: scale !! divisor applied to range, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_stat(x, axis_value, stat_kind_iqr, reduced, &
            reduced_shape, ok, rng=rng, scale=scale, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function iqr_keepdims_rank9_assumed

end module scifort_descriptive_highrank_stats_5_9
