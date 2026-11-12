! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_circular_1_8
    use scifort_descriptive_circular_engine, only : circular_all_assumed_rank, &
        circular_kind_mean, circular_kind_std, circular_kind_var, &
        circular_reduce_assumed_rank, circular_vector
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface circmean
        module procedure circmean_all_rank1
        module procedure circmean_axis_rank1
        module procedure circmean_all_rank2
        module procedure circmean_axis_rank2
        module procedure circmean_all_rank3
        module procedure circmean_axis_rank3
        module procedure circmean_all_rank4
        module procedure circmean_axis_rank4
        module procedure circmean_all_rank5
        module procedure circmean_axis_rank5
        module procedure circmean_all_rank6
        module procedure circmean_axis_rank6
        module procedure circmean_all_rank7
        module procedure circmean_axis_rank7
        module procedure circmean_all_rank8
        module procedure circmean_axis_rank8
    end interface circmean

    interface circmean_all
        module procedure circmean_all_rank1
        module procedure circmean_all_rank2
        module procedure circmean_all_rank3
        module procedure circmean_all_rank4
        module procedure circmean_all_rank5
        module procedure circmean_all_rank6
        module procedure circmean_all_rank7
        module procedure circmean_all_rank8
    end interface circmean_all

    interface circmean_keepdims
        module procedure circmean_keepdims_rank1
        module procedure circmean_keepdims_rank2
        module procedure circmean_keepdims_rank3
        module procedure circmean_keepdims_rank4
        module procedure circmean_keepdims_rank5
        module procedure circmean_keepdims_rank6
        module procedure circmean_keepdims_rank7
        module procedure circmean_keepdims_rank8
    end interface circmean_keepdims

    interface circvar
        module procedure circvar_all_rank1
        module procedure circvar_axis_rank1
        module procedure circvar_all_rank2
        module procedure circvar_axis_rank2
        module procedure circvar_all_rank3
        module procedure circvar_axis_rank3
        module procedure circvar_all_rank4
        module procedure circvar_axis_rank4
        module procedure circvar_all_rank5
        module procedure circvar_axis_rank5
        module procedure circvar_all_rank6
        module procedure circvar_axis_rank6
        module procedure circvar_all_rank7
        module procedure circvar_axis_rank7
        module procedure circvar_all_rank8
        module procedure circvar_axis_rank8
    end interface circvar

    interface circvar_all
        module procedure circvar_all_rank1
        module procedure circvar_all_rank2
        module procedure circvar_all_rank3
        module procedure circvar_all_rank4
        module procedure circvar_all_rank5
        module procedure circvar_all_rank6
        module procedure circvar_all_rank7
        module procedure circvar_all_rank8
    end interface circvar_all

    interface circvar_keepdims
        module procedure circvar_keepdims_rank1
        module procedure circvar_keepdims_rank2
        module procedure circvar_keepdims_rank3
        module procedure circvar_keepdims_rank4
        module procedure circvar_keepdims_rank5
        module procedure circvar_keepdims_rank6
        module procedure circvar_keepdims_rank7
        module procedure circvar_keepdims_rank8
    end interface circvar_keepdims

    interface circstd
        module procedure circstd_all_rank1
        module procedure circstd_axis_rank1
        module procedure circstd_all_rank2
        module procedure circstd_axis_rank2
        module procedure circstd_all_rank3
        module procedure circstd_axis_rank3
        module procedure circstd_all_rank4
        module procedure circstd_axis_rank4
        module procedure circstd_all_rank5
        module procedure circstd_axis_rank5
        module procedure circstd_all_rank6
        module procedure circstd_axis_rank6
        module procedure circstd_all_rank7
        module procedure circstd_axis_rank7
        module procedure circstd_all_rank8
        module procedure circstd_axis_rank8
    end interface circstd

    interface circstd_all
        module procedure circstd_all_rank1
        module procedure circstd_all_rank2
        module procedure circstd_all_rank3
        module procedure circstd_all_rank4
        module procedure circstd_all_rank5
        module procedure circstd_all_rank6
        module procedure circstd_all_rank7
        module procedure circstd_all_rank8
    end interface circstd_all

    interface circstd_keepdims
        module procedure circstd_keepdims_rank1
        module procedure circstd_keepdims_rank2
        module procedure circstd_keepdims_rank3
        module procedure circstd_keepdims_rank4
        module procedure circstd_keepdims_rank5
        module procedure circstd_keepdims_rank6
        module procedure circstd_keepdims_rank7
        module procedure circstd_keepdims_rank8
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

    pure function circmean_all_rank1(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank1

    pure function circmean_axis_rank1(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        if (axis /= 1) then
            values = quiet_nan(0.0_dp)
            return
        end if
        values = circular_vector(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_axis_rank1

    pure function circmean_keepdims_rank1(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        allocate(values(1))
        values(1) = circular_vector(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_keepdims_rank1

    pure function circvar_all_rank1(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank1

    pure function circvar_axis_rank1(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        if (axis /= 1) then
            values = quiet_nan(0.0_dp)
            return
        end if
        values = circular_vector(x, circular_kind_var, high, low, nan_policy)
    end function circvar_axis_rank1

    pure function circvar_keepdims_rank1(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        allocate(values(1))
        values(1) = circular_vector(x, circular_kind_var, high, low, nan_policy)
    end function circvar_keepdims_rank1

    pure function circstd_all_rank1(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank1

    pure function circstd_axis_rank1(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: values

        if (axis /= 1) then
            values = quiet_nan(0.0_dp)
            return
        end if
        values = circular_vector(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_axis_rank1

    pure function circstd_keepdims_rank1(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:)

        integer :: axis_value

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value /= 1) then
            allocate(values(0))
            return
        end if
        allocate(values(1))
        values(1) = circular_vector(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_keepdims_rank1

    pure function circmean_all_rank2(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank2

    pure function circmean_axis_rank2(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank2

    pure function circmean_keepdims_rank2(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank2

    pure function circvar_all_rank2(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank2

    pure function circvar_axis_rank2(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank2

    pure function circvar_keepdims_rank2(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank2

    pure function circstd_all_rank2(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank2

    pure function circstd_axis_rank2(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:)

        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank2

    pure function circstd_keepdims_rank2(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: reduced_shape(1)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank2

    pure function circmean_all_rank3(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank3

    pure function circmean_axis_rank3(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank3

    pure function circmean_keepdims_rank3(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank3

    pure function circvar_all_rank3(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank3

    pure function circvar_axis_rank3(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :)

        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank3

    pure function circvar_keepdims_rank3(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank3

    pure function circstd_all_rank3(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank3

    pure function circstd_axis_rank3(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :)

        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank3

    pure function circstd_keepdims_rank3(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: reduced_shape(2)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank3

    pure function circmean_all_rank4(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank4

    pure function circmean_axis_rank4(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank4

    pure function circmean_keepdims_rank4(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank4

    pure function circvar_all_rank4(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank4

    pure function circvar_axis_rank4(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :)

        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank4

    pure function circvar_keepdims_rank4(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank4

    pure function circstd_all_rank4(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank4

    pure function circstd_axis_rank4(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :)

        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank4

    pure function circstd_keepdims_rank4(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: reduced_shape(3)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank4

    pure function circmean_all_rank5(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank5

    pure function circmean_axis_rank5(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank5

    pure function circmean_keepdims_rank5(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank5

    pure function circvar_all_rank5(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank5

    pure function circvar_axis_rank5(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :)

        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank5

    pure function circvar_keepdims_rank5(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank5

    pure function circstd_all_rank5(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank5

    pure function circstd_axis_rank5(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :)

        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank5

    pure function circstd_keepdims_rank5(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: reduced_shape(4)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank5

    pure function circmean_all_rank6(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank6

    pure function circmean_axis_rank6(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank6

    pure function circmean_keepdims_rank6(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank6

    pure function circvar_all_rank6(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank6

    pure function circvar_axis_rank6(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank6

    pure function circvar_keepdims_rank6(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank6

    pure function circstd_all_rank6(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank6

    pure function circstd_axis_rank6(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank6

    pure function circstd_keepdims_rank6(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: reduced_shape(5)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank6

    pure function circmean_all_rank7(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank7

    pure function circmean_axis_rank7(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank7

    pure function circmean_keepdims_rank7(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank7

    pure function circvar_all_rank7(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank7

    pure function circvar_axis_rank7(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank7

    pure function circvar_keepdims_rank7(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank7

    pure function circstd_all_rank7(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank7

    pure function circstd_axis_rank7(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank7

    pure function circstd_keepdims_rank7(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: reduced_shape(6)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank7

    pure function circmean_all_rank8(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_mean, high, low, nan_policy)
    end function circmean_all_rank8

    pure function circmean_axis_rank8(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_mean, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circmean_axis_rank8

    pure function circmean_keepdims_rank8(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_mean, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circmean_keepdims_rank8

    pure function circvar_all_rank8(x, high, low, nan_policy) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_var, high, low, nan_policy)
    end function circvar_all_rank8

    pure function circvar_axis_rank8(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_var, reduced, &
            out_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circvar_axis_rank8

    pure function circvar_keepdims_rank8(x, high, low, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_var, reduced, &
            reduced_shape, ok, high, low, nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circvar_keepdims_rank8

    pure function circstd_all_rank8(x, high, low, nan_policy, normalize) result(value)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp) :: value

        value = circular_all_assumed_rank(x, circular_kind_std, high, low, nan_policy, normalize)
    end function circstd_all_rank8

    pure function circstd_axis_rank8(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in) :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        call circular_reduce_assumed_rank(x, axis, circular_kind_std, reduced, &
            out_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function circstd_axis_rank8

    pure function circstd_keepdims_rank8(x, high, low, axis, nan_policy, normalize) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! angle observations
        real(dp), intent(in), optional :: high !! upper boundary of a full angle
        real(dp), intent(in), optional :: low !! lower boundary of a full angle
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        logical, intent(in), optional :: normalize !! return unscaled circular std when true
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: reduced_shape(7)
        real(dp), allocatable :: reduced(:)
        logical :: ok

        axis_value = 1
        if (present(axis)) axis_value = axis
        if (axis_value < 1 .or. axis_value > rank(x)) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        call circular_reduce_assumed_rank(x, axis_value, circular_kind_std, reduced, &
            reduced_shape, ok, high, low, nan_policy, normalize)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        full_shape = shape(x)
        full_shape(axis_value) = 1
        values = reshape(reduced, full_shape)
    end function circstd_keepdims_rank8

end module scifort_descriptive_circular_1_8
