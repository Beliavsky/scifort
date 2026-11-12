! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_proportion_trim_9_15
    use scifort_descriptive_proportion_trim_engine, only : &
        proportion_trim_all_assumed_rank, proportion_trim_assumed_rank, &
        trim_mean_all_assumed_rank, trim_mean_reduce_assumed_rank, &
        trim_mean_vector, trim_mode_both, trim_mode_from_tail
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface trim_mean
        module procedure trim_mean_rank9
        module procedure trim_mean_rank10
        module procedure trim_mean_rank11
        module procedure trim_mean_rank12
        module procedure trim_mean_rank13
        module procedure trim_mean_rank14
        module procedure trim_mean_rank15
    end interface trim_mean

    interface trim_mean_all
        module procedure trim_mean_all_rank9
        module procedure trim_mean_all_rank10
        module procedure trim_mean_all_rank11
        module procedure trim_mean_all_rank12
        module procedure trim_mean_all_rank13
        module procedure trim_mean_all_rank14
        module procedure trim_mean_all_rank15
    end interface trim_mean_all

    interface trim_mean_keepdims
        module procedure trim_mean_keepdims_rank9
        module procedure trim_mean_keepdims_rank10
        module procedure trim_mean_keepdims_rank11
        module procedure trim_mean_keepdims_rank12
        module procedure trim_mean_keepdims_rank13
        module procedure trim_mean_keepdims_rank14
        module procedure trim_mean_keepdims_rank15
    end interface trim_mean_keepdims

    interface trimboth
        module procedure trimboth_rank9
        module procedure trimboth_rank10
        module procedure trimboth_rank11
        module procedure trimboth_rank12
        module procedure trimboth_rank13
        module procedure trimboth_rank14
        module procedure trimboth_rank15
    end interface trimboth

    interface trimboth_all
        module procedure trimboth_all_rank9
        module procedure trimboth_all_rank10
        module procedure trimboth_all_rank11
        module procedure trimboth_all_rank12
        module procedure trimboth_all_rank13
        module procedure trimboth_all_rank14
        module procedure trimboth_all_rank15
    end interface trimboth_all

    interface trim1
        module procedure trim1_rank9
        module procedure trim1_rank10
        module procedure trim1_rank11
        module procedure trim1_rank12
        module procedure trim1_rank13
        module procedure trim1_rank14
        module procedure trim1_rank15
    end interface trim1

    interface trim1_all
        module procedure trim1_all_rank9
        module procedure trim1_all_rank10
        module procedure trim1_all_rank11
        module procedure trim1_all_rank12
        module procedure trim1_all_rank13
        module procedure trim1_all_rank14
        module procedure trim1_all_rank15
    end interface trim1_all

    public :: trim_mean
    public :: trim_mean_all
    public :: trim_mean_keepdims
    public :: trimboth
    public :: trimboth_all
    public :: trim1
    public :: trim1_all

contains

    pure function trim_mean_all_rank9(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank9

    pure function trim_mean_rank9(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank9

    pure function trim_mean_keepdims_rank9(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank9

    pure function trimboth_all_rank9(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank9

    pure function trimboth_rank9(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(9)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank9

    pure function trim1_all_rank9(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank9

    pure function trim1_rank9(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(9)
        real(dp), allocatable :: trimmed(:)
        integer :: trim_mode
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank9

    pure function trim_mean_all_rank10(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank10

    pure function trim_mean_rank10(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank10

    pure function trim_mean_keepdims_rank10(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank10

    pure function trimboth_all_rank10(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank10

    pure function trimboth_rank10(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(10)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank10

    pure function trim1_all_rank10(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank10

    pure function trim1_rank10(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(10)
        real(dp), allocatable :: trimmed(:)
        integer :: trim_mode
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank10

    pure function trim_mean_all_rank11(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank11

    pure function trim_mean_rank11(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank11

    pure function trim_mean_keepdims_rank11(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank11

    pure function trimboth_all_rank11(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank11

    pure function trimboth_rank11(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(11)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank11

    pure function trim1_all_rank11(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank11

    pure function trim1_rank11(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(11)
        real(dp), allocatable :: trimmed(:)
        integer :: trim_mode
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank11

    pure function trim_mean_all_rank12(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank12

    pure function trim_mean_rank12(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank12

    pure function trim_mean_keepdims_rank12(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank12

    pure function trimboth_all_rank12(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank12

    pure function trimboth_rank12(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(12)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank12

    pure function trim1_all_rank12(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank12

    pure function trim1_rank12(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(12)
        real(dp), allocatable :: trimmed(:)
        integer :: trim_mode
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank12

    pure function trim_mean_all_rank13(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank13

    pure function trim_mean_rank13(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank13

    pure function trim_mean_keepdims_rank13(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank13

    pure function trimboth_all_rank13(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank13

    pure function trimboth_rank13(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(13)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank13

    pure function trim1_all_rank13(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank13

    pure function trim1_rank13(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(13)
        real(dp), allocatable :: trimmed(:)
        integer :: trim_mode
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank13

    pure function trim_mean_all_rank14(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank14

    pure function trim_mean_rank14(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank14

    pure function trim_mean_keepdims_rank14(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank14

    pure function trimboth_all_rank14(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank14

    pure function trimboth_rank14(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(14)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank14

    pure function trim1_all_rank14(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank14

    pure function trim1_rank14(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(14)
        real(dp), allocatable :: trimmed(:)
        integer :: trim_mode
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank14

    pure function trim_mean_all_rank15(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank15

    pure function trim_mean_rank15(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank15

    pure function trim_mean_keepdims_rank15(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank15

    pure function trimboth_all_rank15(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank15

    pure function trimboth_rank15(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(15)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank15

    pure function trim1_all_rank15(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank15

    pure function trim1_rank15(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(15)
        real(dp), allocatable :: trimmed(:)
        integer :: trim_mode
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank15

end module scifort_descriptive_proportion_trim_9_15
