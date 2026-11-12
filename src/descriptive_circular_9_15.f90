! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_circular_9_15
    use scifort_descriptive_circular_engine, only : circular_all_assumed_rank, &
        circular_kind_mean, circular_kind_std, circular_kind_var, &
        circular_reduce_assumed_rank, circular_vector
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface circmean
        module procedure circmean_all_rank9
        module procedure circmean_axis_rank9
        module procedure circmean_all_rank10
        module procedure circmean_axis_rank10
        module procedure circmean_all_rank11
        module procedure circmean_axis_rank11
        module procedure circmean_all_rank12
        module procedure circmean_axis_rank12
        module procedure circmean_all_rank13
        module procedure circmean_axis_rank13
        module procedure circmean_all_rank14
        module procedure circmean_axis_rank14
        module procedure circmean_all_rank15
        module procedure circmean_axis_rank15
    end interface circmean

    interface circmean_all
        module procedure circmean_all_rank9
        module procedure circmean_all_rank10
        module procedure circmean_all_rank11
        module procedure circmean_all_rank12
        module procedure circmean_all_rank13
        module procedure circmean_all_rank14
        module procedure circmean_all_rank15
    end interface circmean_all

    interface circmean_keepdims
        module procedure circmean_keepdims_rank9
        module procedure circmean_keepdims_rank10
        module procedure circmean_keepdims_rank11
        module procedure circmean_keepdims_rank12
        module procedure circmean_keepdims_rank13
        module procedure circmean_keepdims_rank14
        module procedure circmean_keepdims_rank15
    end interface circmean_keepdims

    interface circvar
        module procedure circvar_all_rank9
        module procedure circvar_axis_rank9
        module procedure circvar_all_rank10
        module procedure circvar_axis_rank10
        module procedure circvar_all_rank11
        module procedure circvar_axis_rank11
        module procedure circvar_all_rank12
        module procedure circvar_axis_rank12
        module procedure circvar_all_rank13
        module procedure circvar_axis_rank13
        module procedure circvar_all_rank14
        module procedure circvar_axis_rank14
        module procedure circvar_all_rank15
        module procedure circvar_axis_rank15
    end interface circvar

    interface circvar_all
        module procedure circvar_all_rank9
        module procedure circvar_all_rank10
        module procedure circvar_all_rank11
        module procedure circvar_all_rank12
        module procedure circvar_all_rank13
        module procedure circvar_all_rank14
        module procedure circvar_all_rank15
    end interface circvar_all

    interface circvar_keepdims
        module procedure circvar_keepdims_rank9
        module procedure circvar_keepdims_rank10
        module procedure circvar_keepdims_rank11
        module procedure circvar_keepdims_rank12
        module procedure circvar_keepdims_rank13
        module procedure circvar_keepdims_rank14
        module procedure circvar_keepdims_rank15
    end interface circvar_keepdims

    interface circstd
        module procedure circstd_all_rank9
        module procedure circstd_axis_rank9
        module procedure circstd_all_rank10
        module procedure circstd_axis_rank10
        module procedure circstd_all_rank11
        module procedure circstd_axis_rank11
        module procedure circstd_all_rank12
        module procedure circstd_axis_rank12
        module procedure circstd_all_rank13
        module procedure circstd_axis_rank13
        module procedure circstd_all_rank14
        module procedure circstd_axis_rank14
        module procedure circstd_all_rank15
        module procedure circstd_axis_rank15
    end interface circstd

    interface circstd_all
        module procedure circstd_all_rank9
        module procedure circstd_all_rank10
        module procedure circstd_all_rank11
        module procedure circstd_all_rank12
        module procedure circstd_all_rank13
        module procedure circstd_all_rank14
        module procedure circstd_all_rank15
    end interface circstd_all

    interface circstd_keepdims
        module procedure circstd_keepdims_rank9
        module procedure circstd_keepdims_rank10
        module procedure circstd_keepdims_rank11
        module procedure circstd_keepdims_rank12
        module procedure circstd_keepdims_rank13
        module procedure circstd_keepdims_rank14
        module procedure circstd_keepdims_rank15
    end interface circstd_keepdims

    public :: circmean
    public :: circmean_all
    public :: circmean_keepdims
    public :: circvar
    public :: circvar_all
    public :: circvar_keepdims
    public :: circstd
    public :: circstd_all
    public :: circstd_keepdims

contains

    pure function circmean_all_rank9(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank9

    pure function circmean_axis_rank9(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank9

    pure function circmean_keepdims_rank9(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank9

    pure function circvar_all_rank9(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank9

    pure function circvar_axis_rank9(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank9

    pure function circvar_keepdims_rank9(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank9

    pure function circstd_all_rank9(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank9

    pure function circstd_axis_rank9(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank9

    pure function circstd_keepdims_rank9(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: reduced_shape(8)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank9

    pure function circmean_all_rank10(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank10

    pure function circmean_axis_rank10(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank10

    pure function circmean_keepdims_rank10(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank10

    pure function circvar_all_rank10(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank10

    pure function circvar_axis_rank10(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank10

    pure function circvar_keepdims_rank10(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank10

    pure function circstd_all_rank10(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank10

    pure function circstd_axis_rank10(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank10

    pure function circstd_keepdims_rank10(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: reduced_shape(9)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank10

    pure function circmean_all_rank11(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank11

    pure function circmean_axis_rank11(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank11

    pure function circmean_keepdims_rank11(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank11

    pure function circvar_all_rank11(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank11

    pure function circvar_axis_rank11(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank11

    pure function circvar_keepdims_rank11(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank11

    pure function circstd_all_rank11(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank11

    pure function circstd_axis_rank11(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank11

    pure function circstd_keepdims_rank11(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: reduced_shape(10)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank11

    pure function circmean_all_rank12(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank12

    pure function circmean_axis_rank12(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank12

    pure function circmean_keepdims_rank12(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank12

    pure function circvar_all_rank12(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank12

    pure function circvar_axis_rank12(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank12

    pure function circvar_keepdims_rank12(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank12

    pure function circstd_all_rank12(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank12

    pure function circstd_axis_rank12(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank12

    pure function circstd_keepdims_rank12(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: reduced_shape(11)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank12

    pure function circmean_all_rank13(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank13

    pure function circmean_axis_rank13(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank13

    pure function circmean_keepdims_rank13(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank13

    pure function circvar_all_rank13(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank13

    pure function circvar_axis_rank13(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank13

    pure function circvar_keepdims_rank13(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank13

    pure function circstd_all_rank13(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank13

    pure function circstd_axis_rank13(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank13

    pure function circstd_keepdims_rank13(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: reduced_shape(12)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank13

    pure function circmean_all_rank14(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank14

    pure function circmean_axis_rank14(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank14

    pure function circmean_keepdims_rank14(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank14

    pure function circvar_all_rank14(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank14

    pure function circvar_axis_rank14(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank14

    pure function circvar_keepdims_rank14(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank14

    pure function circstd_all_rank14(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank14

    pure function circstd_axis_rank14(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank14

    pure function circstd_keepdims_rank14(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: reduced_shape(13)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank14

    pure function circmean_all_rank15(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank15

    pure function circmean_axis_rank15(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank15

    pure function circmean_keepdims_rank15(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank15

    pure function circvar_all_rank15(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank15

    pure function circvar_axis_rank15(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank15

    pure function circvar_keepdims_rank15(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank15

    pure function circstd_all_rank15(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank15

    pure function circstd_axis_rank15(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank15

    pure function circstd_keepdims_rank15(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: reduced_shape(14)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank15

end module scifort_descriptive_circular_9_15
