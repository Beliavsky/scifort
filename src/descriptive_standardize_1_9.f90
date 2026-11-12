! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_standardize_1_9
    use scifort_descriptive_standardize_engine, only : gstd_all_assumed_rank, &
        gstd_reduce_assumed_rank, gstd_reduce_dims_keepdims_assumed_rank, &
        gzscore_all_assumed_rank, gzscore_assumed_rank, &
        zmap_all_assumed_rank, zmap_assumed_rank, zscore_all_assumed_rank, &
        zscore_assumed_rank
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface zscore_all
        module procedure zscore_all_rank1
        module procedure zscore_all_rank2
        module procedure zscore_all_rank3
        module procedure zscore_all_rank4
        module procedure zscore_all_rank5
        module procedure zscore_all_rank6
        module procedure zscore_all_rank7
        module procedure zscore_all_rank8
        module procedure zscore_all_rank9
    end interface zscore_all

    interface zmap_all
        module procedure zmap_all_rank1
        module procedure zmap_all_rank2
        module procedure zmap_all_rank3
        module procedure zmap_all_rank4
        module procedure zmap_all_rank5
        module procedure zmap_all_rank6
        module procedure zmap_all_rank7
        module procedure zmap_all_rank8
        module procedure zmap_all_rank9
    end interface zmap_all

    interface gzscore_all
        module procedure gzscore_all_rank1
        module procedure gzscore_all_rank2
        module procedure gzscore_all_rank3
        module procedure gzscore_all_rank4
        module procedure gzscore_all_rank5
        module procedure gzscore_all_rank6
        module procedure gzscore_all_rank7
        module procedure gzscore_all_rank8
        module procedure gzscore_all_rank9
    end interface gzscore_all

    interface gstd_all
        module procedure gstd_all_rank1
        module procedure gstd_all_rank2
        module procedure gstd_all_rank3
        module procedure gstd_all_rank4
        module procedure gstd_all_rank5
        module procedure gstd_all_rank6
        module procedure gstd_all_rank7
        module procedure gstd_all_rank8
        module procedure gstd_all_rank9
    end interface gstd_all

    interface zscore
        module procedure zscore_rank1
        module procedure zscore_rank2
        module procedure zscore_rank3
        module procedure zscore_rank4
        module procedure zscore_rank5
        module procedure zscore_rank6
        module procedure zscore_rank7
        module procedure zscore_rank8
        module procedure zscore_rank9
    end interface zscore

    interface zmap
        module procedure zmap_rank1
        module procedure zmap_rank2
        module procedure zmap_rank3
        module procedure zmap_rank4
        module procedure zmap_rank5
        module procedure zmap_rank6
        module procedure zmap_rank7
        module procedure zmap_rank8
        module procedure zmap_rank9
    end interface zmap

    interface gzscore
        module procedure gzscore_rank1
        module procedure gzscore_rank2
        module procedure gzscore_rank3
        module procedure gzscore_rank4
        module procedure gzscore_rank5
        module procedure gzscore_rank6
        module procedure gzscore_rank7
        module procedure gzscore_rank8
        module procedure gzscore_rank9
    end interface gzscore

    interface gstd
        module procedure gstd_rank1
        module procedure gstd_rank2
        module procedure gstd_rank3
        module procedure gstd_rank4
        module procedure gstd_rank5
        module procedure gstd_rank6
        module procedure gstd_rank7
        module procedure gstd_rank8
        module procedure gstd_rank9
    end interface gstd

    interface gstd_keepdims
        module procedure gstd_keepdims_rank1
        module procedure gstd_keepdims_rank2
        module procedure gstd_keepdims_rank3
        module procedure gstd_keepdims_rank4
        module procedure gstd_keepdims_rank5
        module procedure gstd_keepdims_rank6
        module procedure gstd_keepdims_rank7
        module procedure gstd_keepdims_rank8
        module procedure gstd_keepdims_rank9
        module procedure gstd_keepdims_dims_rank1
        module procedure gstd_keepdims_dims_rank2
        module procedure gstd_keepdims_dims_rank3
        module procedure gstd_keepdims_dims_rank4
        module procedure gstd_keepdims_dims_rank5
        module procedure gstd_keepdims_dims_rank6
        module procedure gstd_keepdims_dims_rank7
        module procedure gstd_keepdims_dims_rank8
        module procedure gstd_keepdims_dims_rank9
    end interface gstd_keepdims

    public :: zscore_all
    public :: zmap_all
    public :: gzscore_all
    public :: gstd_all
    public :: zscore
    public :: zmap
    public :: gzscore
    public :: gstd
    public :: gstd_keepdims

contains

    pure function zmap_rank1(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:) !! values standardized against compare
        real(dp), intent(in) :: compare(:) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 1; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(1)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank1

    pure function zmap_rank2(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(2)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank2

    pure function zmap_rank3(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(3)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank3

    pure function zmap_rank4(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(4)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank4

    pure function zmap_rank5(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(5)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank5

    pure function zmap_rank6(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(6)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank6

    pure function zmap_rank7(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(7)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank7

    pure function zmap_rank8(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(8)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank8

    pure function zmap_rank9(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(9)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank9

    pure function zscore_rank1(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 1; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank1

    pure function gzscore_rank1(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 1; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank1

    pure function gstd_rank1(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 1; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        integer :: axis_value
        integer :: reduced_shape(0)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            value = quiet_nan(0.0_dp)
            return
        end if
        value = reduced(1)
    end function gstd_rank1

    pure function gstd_keepdims_rank1(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 1; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        integer :: reduced_shape(0)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank1

    pure function zscore_rank2(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank2

    pure function gzscore_rank2(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank2

    pure function gstd_rank2(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:)

        integer :: axis_value
        integer :: reduced_shape(1)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank2

    pure function gstd_keepdims_rank2(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 2; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank2

    pure function zscore_rank3(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank3

    pure function gzscore_rank3(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank3

    pure function gstd_rank3(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :)

        integer :: axis_value
        integer :: reduced_shape(2)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank3

    pure function gstd_keepdims_rank3(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 3; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank3

    pure function zscore_rank4(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank4

    pure function gzscore_rank4(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank4

    pure function gstd_rank4(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :)

        integer :: axis_value
        integer :: reduced_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank4

    pure function gstd_keepdims_rank4(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank4

    pure function zscore_rank5(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank5

    pure function gzscore_rank5(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank5

    pure function gstd_rank5(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank5

    pure function gstd_keepdims_rank5(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank5

    pure function zscore_rank6(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank6

    pure function gzscore_rank6(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank6

    pure function gstd_rank6(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank6

    pure function gstd_keepdims_rank6(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank6

    pure function zscore_rank7(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank7

    pure function gzscore_rank7(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank7

    pure function gstd_rank7(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank7

    pure function gstd_keepdims_rank7(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank7

    pure function zscore_rank8(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank8

    pure function gzscore_rank8(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank8

    pure function gstd_rank8(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank8

    pure function gstd_keepdims_rank8(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank8

    pure function zscore_rank9(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank9

    pure function gzscore_rank9(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank9

    pure function gstd_rank9(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank9

    pure function gstd_keepdims_rank9(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank9

    pure function zmap_all_rank1(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(1)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank1
    pure function zscore_all_rank1(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(1)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank1
    pure function gzscore_all_rank1(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(1)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank1
    pure function gstd_all_rank1(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank1
    pure function zmap_all_rank2(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(2)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank2
    pure function zscore_all_rank2(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(2)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank2
    pure function gzscore_all_rank2(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(2)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank2
    pure function gstd_all_rank2(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank2
    pure function zmap_all_rank3(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(3)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank3
    pure function zscore_all_rank3(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(3)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank3
    pure function gzscore_all_rank3(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(3)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank3
    pure function gstd_all_rank3(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank3
    pure function zmap_all_rank4(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(4)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank4
    pure function zscore_all_rank4(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(4)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank4
    pure function gzscore_all_rank4(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(4)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank4
    pure function gstd_all_rank4(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank4
    pure function zmap_all_rank5(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(5)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank5
    pure function zscore_all_rank5(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(5)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank5
    pure function gzscore_all_rank5(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(5)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank5
    pure function gstd_all_rank5(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank5
    pure function zmap_all_rank6(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(6)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank6
    pure function zscore_all_rank6(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(6)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank6
    pure function gzscore_all_rank6(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(6)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank6
    pure function gstd_all_rank6(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank6
    pure function zmap_all_rank7(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(7)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank7
    pure function zscore_all_rank7(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(7)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank7
    pure function gzscore_all_rank7(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(7)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank7
    pure function gstd_all_rank7(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank7
    pure function zmap_all_rank8(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, &
            :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, &
            :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(8)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank8
    pure function zscore_all_rank8(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(8)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank8
    pure function gzscore_all_rank8(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(8)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank8
    pure function gstd_all_rank8(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank8
    pure function zmap_all_rank9(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, &
            :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, &
            :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(9)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank9
    pure function zscore_all_rank9(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(9)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank9
    pure function gzscore_all_rank9(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(9)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank9
    pure function gstd_all_rank9(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank9

    pure function gstd_keepdims_dims_rank1(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:)

        integer :: full_shape(1)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank1
    pure function gstd_keepdims_dims_rank2(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :)

        integer :: full_shape(2)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank2
    pure function gstd_keepdims_dims_rank3(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: full_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank3
    pure function gstd_keepdims_dims_rank4(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: full_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank4
    pure function gstd_keepdims_dims_rank5(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: full_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank5
    pure function gstd_keepdims_dims_rank6(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: full_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank6
    pure function gstd_keepdims_dims_rank7(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: full_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank7
    pure function gstd_keepdims_dims_rank8(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: full_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank8
    pure function gstd_keepdims_dims_rank9(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: full_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank9

end module scifort_descriptive_standardize_1_9
