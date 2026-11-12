! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_standardize_10_15
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
        module procedure zscore_all_rank10
        module procedure zscore_all_rank11
        module procedure zscore_all_rank12
        module procedure zscore_all_rank13
        module procedure zscore_all_rank14
        module procedure zscore_all_rank15
    end interface zscore_all

    interface zmap_all
        module procedure zmap_all_rank10
        module procedure zmap_all_rank11
        module procedure zmap_all_rank12
        module procedure zmap_all_rank13
        module procedure zmap_all_rank14
        module procedure zmap_all_rank15
    end interface zmap_all

    interface gzscore_all
        module procedure gzscore_all_rank10
        module procedure gzscore_all_rank11
        module procedure gzscore_all_rank12
        module procedure gzscore_all_rank13
        module procedure gzscore_all_rank14
        module procedure gzscore_all_rank15
    end interface gzscore_all

    interface gstd_all
        module procedure gstd_all_rank10
        module procedure gstd_all_rank11
        module procedure gstd_all_rank12
        module procedure gstd_all_rank13
        module procedure gstd_all_rank14
        module procedure gstd_all_rank15
    end interface gstd_all

    interface zscore
        module procedure zscore_rank10
        module procedure zscore_rank11
        module procedure zscore_rank12
        module procedure zscore_rank13
        module procedure zscore_rank14
        module procedure zscore_rank15
    end interface zscore

    interface zmap
        module procedure zmap_rank10
        module procedure zmap_rank11
        module procedure zmap_rank12
        module procedure zmap_rank13
        module procedure zmap_rank14
        module procedure zmap_rank15
    end interface zmap

    interface gzscore
        module procedure gzscore_rank10
        module procedure gzscore_rank11
        module procedure gzscore_rank12
        module procedure gzscore_rank13
        module procedure gzscore_rank14
        module procedure gzscore_rank15
    end interface gzscore

    interface gstd
        module procedure gstd_rank10
        module procedure gstd_rank11
        module procedure gstd_rank12
        module procedure gstd_rank13
        module procedure gstd_rank14
        module procedure gstd_rank15
    end interface gstd

    interface gstd_keepdims
        module procedure gstd_keepdims_rank10
        module procedure gstd_keepdims_rank11
        module procedure gstd_keepdims_rank12
        module procedure gstd_keepdims_rank13
        module procedure gstd_keepdims_rank14
        module procedure gstd_keepdims_rank15
        module procedure gstd_keepdims_dims_rank10
        module procedure gstd_keepdims_dims_rank11
        module procedure gstd_keepdims_dims_rank12
        module procedure gstd_keepdims_dims_rank13
        module procedure gstd_keepdims_dims_rank14
        module procedure gstd_keepdims_dims_rank15
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

    pure function zmap_rank10(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(10)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank10

    pure function zmap_rank11(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(11)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank11

    pure function zmap_rank12(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(12)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank12

    pure function zmap_rank13(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(13)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank13

    pure function zmap_rank14(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(14)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank14

    pure function zmap_rank15(scores, compare, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! values standardized against compare
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        integer :: final_shape(15)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zmap_assumed_rank(scores, compare, axis_value, ddof, nan_policy, &
            flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_rank15

    pure function zscore_rank10(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank10

    pure function gzscore_rank10(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank10

    pure function gstd_rank10(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank10

    pure function gstd_keepdims_rank10(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank10

    pure function zscore_rank11(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank11

    pure function gzscore_rank11(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank11

    pure function gstd_rank11(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank11

    pure function gstd_keepdims_rank11(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank11

    pure function zscore_rank12(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank12

    pure function gzscore_rank12(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank12

    pure function gstd_rank12(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank12

    pure function gstd_keepdims_rank12(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank12

    pure function zscore_rank13(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank13

    pure function gzscore_rank13(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank13

    pure function gstd_rank13(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank13

    pure function gstd_keepdims_rank13(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank13

    pure function zscore_rank14(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank14

    pure function gzscore_rank14(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank14

    pure function gstd_rank14(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank14

    pure function gstd_keepdims_rank14(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank14

    pure function zscore_rank15(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call zscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function zscore_rank15

    pure function gzscore_rank15(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer, allocatable :: actual_shape(:)
        real(dp), allocatable :: flat(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gzscore_assumed_rank(x, axis_value, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(flat, shape(x))
    end function gzscore_rank15

    pure function gstd_rank15(x, axis, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: value(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(value(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        value = reshape(reduced, reduced_shape)
    end function gstd_rank15

    pure function gstd_keepdims_rank15(x, axis, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations; standardized along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call gstd_reduce_assumed_rank(x, axis_value, ddof, nan_policy, reduced, &
            reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_rank15

    pure function zmap_all_rank10(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, &
            :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, &
            :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(10)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank10
    pure function zscore_all_rank10(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(10)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank10
    pure function gzscore_all_rank10(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(10)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank10
    pure function gstd_all_rank10(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank10
    pure function zmap_all_rank11(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, &
            :, :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, &
            :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(11)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank11
    pure function zscore_all_rank11(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(11)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank11
    pure function gzscore_all_rank11(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(11)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank11
    pure function gstd_all_rank11(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank11
    pure function zmap_all_rank12(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, &
            :, :, :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, &
            :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(12)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank12
    pure function zscore_all_rank12(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(12)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank12
    pure function gzscore_all_rank12(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(12)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank12
    pure function gstd_all_rank12(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank12
    pure function zmap_all_rank13(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, &
            :, :, :, :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, &
            :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(13)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank13
    pure function zscore_all_rank13(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(13)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank13
    pure function gzscore_all_rank13(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(13)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank13
    pure function gstd_all_rank13(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank13
    pure function zmap_all_rank14(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(14)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank14
    pure function zscore_all_rank14(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(14)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank14
    pure function gzscore_all_rank14(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(14)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank14
    pure function gstd_all_rank14(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank14
    pure function zmap_all_rank15(scores, compare, ddof, nan_policy) result(values)
        real(dp), intent(in) :: scores(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :, &
            :) !! values standardized against whole compare array
        real(dp), intent(in) :: compare(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :, &
            :) !! reference values supplying mean and variance
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(15)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zmap_all_assumed_rank(scores, compare, ddof, nan_policy, flat, &
            actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zmap_all_rank15
    pure function zscore_all_rank15(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :, &
            :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(15)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call zscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function zscore_all_rank15
    pure function gzscore_all_rank15(x, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :, &
            :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 0
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer, allocatable :: actual_shape(:)
        integer :: final_shape(15)
        real(dp), allocatable :: flat(:)
        logical :: ok

        call gzscore_all_assumed_rank(x, ddof, nan_policy, flat, actual_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        final_shape = actual_shape
        values = reshape(flat, final_shape)
    end function gzscore_all_rank15
    pure function gstd_all_rank15(x, ddof, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, &
            :, :, :, :, :, :, :, &
            :) !! observations reduced over every dimension
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp) :: value

        logical :: ok

        call gstd_all_assumed_rank(x, ddof, nan_policy, value, ok)
        if (.not. ok) value = quiet_nan(0.0_dp)
    end function gstd_all_rank15

    pure function gstd_keepdims_dims_rank10(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: full_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank10
    pure function gstd_keepdims_dims_rank11(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: full_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank11
    pure function gstd_keepdims_dims_rank12(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: full_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank12
    pure function gstd_keepdims_dims_rank13(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: full_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank13
    pure function gstd_keepdims_dims_rank14(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: full_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank14
    pure function gstd_keepdims_dims_rank15(x, axes, ddof, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations reduced jointly
        integer, intent(in) :: axes(:) !! unique one-based dimensions reduced jointly
        integer, intent(in), optional :: ddof !! divisor adjustment, default 1
        character(len=*), intent(in), optional :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: full_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        call gstd_reduce_dims_keepdims_assumed_rank(x, axes, ddof, nan_policy, reduced, &
            full_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, full_shape)
    end function gstd_keepdims_dims_rank15

end module scifort_descriptive_standardize_10_15
