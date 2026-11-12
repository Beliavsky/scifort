! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_proportion_trim_1_8
    use scifort_descriptive_proportion_trim_engine, only : &
        proportion_trim_all_assumed_rank, proportion_trim_assumed_rank, &
        trim_mean_all_assumed_rank, trim_mean_reduce_assumed_rank, &
        trim_mean_vector, trim_mode_both, trim_mode_from_tail
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface trim_mean
        module procedure trim_mean_rank1
        module procedure trim_mean_rank2
        module procedure trim_mean_rank3
        module procedure trim_mean_rank4
        module procedure trim_mean_rank5
        module procedure trim_mean_rank6
        module procedure trim_mean_rank7
        module procedure trim_mean_rank8
    end interface trim_mean

    interface trim_mean_all
        module procedure trim_mean_all_rank1
        module procedure trim_mean_all_rank2
        module procedure trim_mean_all_rank3
        module procedure trim_mean_all_rank4
        module procedure trim_mean_all_rank5
        module procedure trim_mean_all_rank6
        module procedure trim_mean_all_rank7
        module procedure trim_mean_all_rank8
    end interface trim_mean_all

    interface trim_mean_keepdims
        module procedure trim_mean_keepdims_rank1
        module procedure trim_mean_keepdims_rank2
        module procedure trim_mean_keepdims_rank3
        module procedure trim_mean_keepdims_rank4
        module procedure trim_mean_keepdims_rank5
        module procedure trim_mean_keepdims_rank6
        module procedure trim_mean_keepdims_rank7
        module procedure trim_mean_keepdims_rank8
    end interface trim_mean_keepdims

    interface trimboth
        module procedure trimboth_rank1
        module procedure trimboth_rank2
        module procedure trimboth_rank3
        module procedure trimboth_rank4
        module procedure trimboth_rank5
        module procedure trimboth_rank6
        module procedure trimboth_rank7
        module procedure trimboth_rank8
    end interface trimboth

    interface trimboth_all
        module procedure trimboth_all_rank1
        module procedure trimboth_all_rank2
        module procedure trimboth_all_rank3
        module procedure trimboth_all_rank4
        module procedure trimboth_all_rank5
        module procedure trimboth_all_rank6
        module procedure trimboth_all_rank7
        module procedure trimboth_all_rank8
    end interface trimboth_all

    interface trim1
        module procedure trim1_rank1
        module procedure trim1_rank2
        module procedure trim1_rank3
        module procedure trim1_rank4
        module procedure trim1_rank5
        module procedure trim1_rank6
        module procedure trim1_rank7
        module procedure trim1_rank8
    end interface trim1

    interface trim1_all
        module procedure trim1_all_rank1
        module procedure trim1_all_rank2
        module procedure trim1_all_rank3
        module procedure trim1_all_rank4
        module procedure trim1_all_rank5
        module procedure trim1_all_rank6
        module procedure trim1_all_rank7
        module procedure trim1_all_rank8
    end interface trim1_all

    public :: trim_mean
    public :: trim_mean_all
    public :: trim_mean_keepdims
    public :: trimboth
    public :: trimboth_all
    public :: trim1
    public :: trim1_all

contains

    pure function trim_mean_all_rank1(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank1

    pure function trim_mean_rank1(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = trim_mean_vector(x, proportiontocut, nan_policy)
    end function trim_mean_rank1

    pure function trim_mean_keepdims_rank1(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
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
        scalar_value = trim_mean_vector(x, proportiontocut, nan_policy)
        values = reshape([scalar_value], full_shape)
    end function trim_mean_keepdims_rank1

    pure function trimboth_all_rank1(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank1

    pure function trimboth_rank1(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: out_shape(1)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank1

    pure function trim1_all_rank1(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank1

    pure function trim1_rank1(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: out_shape(1)
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
            allocate(values(0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank1

    pure function trim_mean_all_rank2(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank2

    pure function trim_mean_rank2(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:)

        integer :: axis_value
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank2

    pure function trim_mean_keepdims_rank2(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank2

    pure function trimboth_all_rank2(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank2

    pure function trimboth_rank2(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: out_shape(2)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank2

    pure function trim1_all_rank2(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank2

    pure function trim1_rank2(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: out_shape(2)
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
            allocate(values(0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank2

    pure function trim_mean_all_rank3(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank3

    pure function trim_mean_rank3(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :)

        integer :: axis_value
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank3

    pure function trim_mean_keepdims_rank3(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank3

    pure function trimboth_all_rank3(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank3

    pure function trimboth_rank3(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: out_shape(3)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank3

    pure function trim1_all_rank3(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank3

    pure function trim1_rank3(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: out_shape(3)
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
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank3

    pure function trim_mean_all_rank4(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank4

    pure function trim_mean_rank4(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :)

        integer :: axis_value
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank4

    pure function trim_mean_keepdims_rank4(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank4

    pure function trimboth_all_rank4(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank4

    pure function trimboth_rank4(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: out_shape(4)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank4

    pure function trim1_all_rank4(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank4

    pure function trim1_rank4(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: out_shape(4)
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
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank4

    pure function trim_mean_all_rank5(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank5

    pure function trim_mean_rank5(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank5

    pure function trim_mean_keepdims_rank5(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank5

    pure function trimboth_all_rank5(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank5

    pure function trimboth_rank5(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(5)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank5

    pure function trim1_all_rank5(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank5

    pure function trim1_rank5(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(5)
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
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank5

    pure function trim_mean_all_rank6(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank6

    pure function trim_mean_rank6(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank6

    pure function trim_mean_keepdims_rank6(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank6

    pure function trimboth_all_rank6(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank6

    pure function trimboth_rank6(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(6)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank6

    pure function trim1_all_rank6(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank6

    pure function trim1_rank6(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(6)
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
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank6

    pure function trim_mean_all_rank7(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank7

    pure function trim_mean_rank7(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank7

    pure function trim_mean_keepdims_rank7(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank7

    pure function trimboth_all_rank7(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank7

    pure function trimboth_rank7(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(7)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank7

    pure function trim1_all_rank7(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank7

    pure function trim1_rank7(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(7)
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
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank7

    pure function trim_mean_all_rank8(x, proportiontocut, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        value = trim_mean_all_assumed_rank(x, proportiontocut, nan_policy)
    end function trim_mean_all_rank8

    pure function trim_mean_rank8(x, proportiontocut, axis, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function trim_mean_rank8

    pure function trim_mean_keepdims_rank8(x, proportiontocut, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call trim_mean_reduce_assumed_rank(x, axis_value, proportiontocut, &
            reduced, reduced_shape, ok, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function trim_mean_keepdims_rank8

    pure function trimboth_all_rank8(x, proportiontocut) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        real(dp), allocatable :: values(:)

        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode_both)
    end function trimboth_all_rank8

    pure function trimboth_rank8(x, proportiontocut, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from each tail
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(8)
        real(dp), allocatable :: trimmed(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call proportion_trim_assumed_rank(x, axis_value, proportiontocut, &
            trim_mode_both, trimmed, out_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trimboth_rank8

    pure function trim1_all_rank8(x, proportiontocut, tail) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        real(dp), allocatable :: values(:)

        integer :: trim_mode

        trim_mode = trim_mode_from_tail('right')
        if (present(tail)) trim_mode = trim_mode_from_tail(tail)
        values = proportion_trim_all_assumed_rank(x, proportiontocut, trim_mode)
    end function trim1_all_rank8

    pure function trim1_rank8(x, proportiontocut, tail, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: proportiontocut !! fraction removed from one tail
        character(len=*), intent(in), optional :: tail !! 'left' or 'right', default right
        integer, intent(in), optional :: axis !! one-based Fortran trimming axis
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: out_shape(8)
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
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(trimmed, out_shape)
    end function trim1_rank8

end module scifort_descriptive_proportion_trim_1_8
