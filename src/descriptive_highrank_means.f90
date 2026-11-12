! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_highrank_means
    use scifort_descriptive_lowrank, only : gmean_low => gmean, hmean_low => hmean, &
        pmean_low => pmean
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter :: mean_kind_gmean = 1
    integer, parameter :: mean_kind_hmean = 2
    integer, parameter :: mean_kind_pmean = 3

    interface gmean
        module procedure gmean_rank5_assumed
        module procedure gmean_rank5_assumed_nan_policy
        module procedure gmean_rank6_assumed
        module procedure gmean_rank6_assumed_nan_policy
        module procedure gmean_rank7_assumed
        module procedure gmean_rank7_assumed_nan_policy
        module procedure gmean_rank8_assumed
        module procedure gmean_rank8_assumed_nan_policy
        module procedure gmean_rank9_assumed
        module procedure gmean_rank9_assumed_nan_policy
        module procedure gmean_rank10_assumed
        module procedure gmean_rank10_assumed_nan_policy
        module procedure gmean_rank11_assumed
        module procedure gmean_rank11_assumed_nan_policy
        module procedure gmean_rank12_assumed
        module procedure gmean_rank12_assumed_nan_policy
        module procedure gmean_rank13_assumed
        module procedure gmean_rank13_assumed_nan_policy
        module procedure gmean_rank14_assumed
        module procedure gmean_rank14_assumed_nan_policy
        module procedure gmean_rank15_assumed
        module procedure gmean_rank15_assumed_nan_policy
    end interface gmean

    interface hmean
        module procedure hmean_rank5_assumed
        module procedure hmean_rank5_assumed_nan_policy
        module procedure hmean_rank6_assumed
        module procedure hmean_rank6_assumed_nan_policy
        module procedure hmean_rank7_assumed
        module procedure hmean_rank7_assumed_nan_policy
        module procedure hmean_rank8_assumed
        module procedure hmean_rank8_assumed_nan_policy
        module procedure hmean_rank9_assumed
        module procedure hmean_rank9_assumed_nan_policy
        module procedure hmean_rank10_assumed
        module procedure hmean_rank10_assumed_nan_policy
        module procedure hmean_rank11_assumed
        module procedure hmean_rank11_assumed_nan_policy
        module procedure hmean_rank12_assumed
        module procedure hmean_rank12_assumed_nan_policy
        module procedure hmean_rank13_assumed
        module procedure hmean_rank13_assumed_nan_policy
        module procedure hmean_rank14_assumed
        module procedure hmean_rank14_assumed_nan_policy
        module procedure hmean_rank15_assumed
        module procedure hmean_rank15_assumed_nan_policy
    end interface hmean

    interface pmean
        module procedure pmean_rank5_assumed
        module procedure pmean_rank5_assumed_nan_policy
        module procedure pmean_rank6_assumed
        module procedure pmean_rank6_assumed_nan_policy
        module procedure pmean_rank7_assumed
        module procedure pmean_rank7_assumed_nan_policy
        module procedure pmean_rank8_assumed
        module procedure pmean_rank8_assumed_nan_policy
        module procedure pmean_rank9_assumed
        module procedure pmean_rank9_assumed_nan_policy
        module procedure pmean_rank10_assumed
        module procedure pmean_rank10_assumed_nan_policy
        module procedure pmean_rank11_assumed
        module procedure pmean_rank11_assumed_nan_policy
        module procedure pmean_rank12_assumed
        module procedure pmean_rank12_assumed_nan_policy
        module procedure pmean_rank13_assumed
        module procedure pmean_rank13_assumed_nan_policy
        module procedure pmean_rank14_assumed
        module procedure pmean_rank14_assumed_nan_policy
        module procedure pmean_rank15_assumed
        module procedure pmean_rank15_assumed_nan_policy
    end interface pmean

    interface gmean_keepdims
        module procedure gmean_rank5_assumed_keepdims
        module procedure gmean_rank5_assumed_keepdims_nan_policy
        module procedure gmean_rank6_assumed_keepdims
        module procedure gmean_rank6_assumed_keepdims_nan_policy
        module procedure gmean_rank7_assumed_keepdims
        module procedure gmean_rank7_assumed_keepdims_nan_policy
        module procedure gmean_rank8_assumed_keepdims
        module procedure gmean_rank8_assumed_keepdims_nan_policy
        module procedure gmean_rank9_assumed_keepdims
        module procedure gmean_rank9_assumed_keepdims_nan_policy
        module procedure gmean_rank10_assumed_keepdims
        module procedure gmean_rank10_assumed_keepdims_nan_policy
        module procedure gmean_rank11_assumed_keepdims
        module procedure gmean_rank11_assumed_keepdims_nan_policy
        module procedure gmean_rank12_assumed_keepdims
        module procedure gmean_rank12_assumed_keepdims_nan_policy
        module procedure gmean_rank13_assumed_keepdims
        module procedure gmean_rank13_assumed_keepdims_nan_policy
        module procedure gmean_rank14_assumed_keepdims
        module procedure gmean_rank14_assumed_keepdims_nan_policy
        module procedure gmean_rank15_assumed_keepdims
        module procedure gmean_rank15_assumed_keepdims_nan_policy
    end interface gmean_keepdims

    interface hmean_keepdims
        module procedure hmean_rank5_assumed_keepdims
        module procedure hmean_rank5_assumed_keepdims_nan_policy
        module procedure hmean_rank6_assumed_keepdims
        module procedure hmean_rank6_assumed_keepdims_nan_policy
        module procedure hmean_rank7_assumed_keepdims
        module procedure hmean_rank7_assumed_keepdims_nan_policy
        module procedure hmean_rank8_assumed_keepdims
        module procedure hmean_rank8_assumed_keepdims_nan_policy
        module procedure hmean_rank9_assumed_keepdims
        module procedure hmean_rank9_assumed_keepdims_nan_policy
        module procedure hmean_rank10_assumed_keepdims
        module procedure hmean_rank10_assumed_keepdims_nan_policy
        module procedure hmean_rank11_assumed_keepdims
        module procedure hmean_rank11_assumed_keepdims_nan_policy
        module procedure hmean_rank12_assumed_keepdims
        module procedure hmean_rank12_assumed_keepdims_nan_policy
        module procedure hmean_rank13_assumed_keepdims
        module procedure hmean_rank13_assumed_keepdims_nan_policy
        module procedure hmean_rank14_assumed_keepdims
        module procedure hmean_rank14_assumed_keepdims_nan_policy
        module procedure hmean_rank15_assumed_keepdims
        module procedure hmean_rank15_assumed_keepdims_nan_policy
    end interface hmean_keepdims

    interface pmean_keepdims
        module procedure pmean_rank5_assumed_keepdims
        module procedure pmean_rank5_assumed_keepdims_nan_policy
        module procedure pmean_rank6_assumed_keepdims
        module procedure pmean_rank6_assumed_keepdims_nan_policy
        module procedure pmean_rank7_assumed_keepdims
        module procedure pmean_rank7_assumed_keepdims_nan_policy
        module procedure pmean_rank8_assumed_keepdims
        module procedure pmean_rank8_assumed_keepdims_nan_policy
        module procedure pmean_rank9_assumed_keepdims
        module procedure pmean_rank9_assumed_keepdims_nan_policy
        module procedure pmean_rank10_assumed_keepdims
        module procedure pmean_rank10_assumed_keepdims_nan_policy
        module procedure pmean_rank11_assumed_keepdims
        module procedure pmean_rank11_assumed_keepdims_nan_policy
        module procedure pmean_rank12_assumed_keepdims
        module procedure pmean_rank12_assumed_keepdims_nan_policy
        module procedure pmean_rank13_assumed_keepdims
        module procedure pmean_rank13_assumed_keepdims_nan_policy
        module procedure pmean_rank14_assumed_keepdims
        module procedure pmean_rank14_assumed_keepdims_nan_policy
        module procedure pmean_rank15_assumed_keepdims
        module procedure pmean_rank15_assumed_keepdims_nan_policy
    end interface pmean_keepdims

    public :: gmean
    public :: gmean_keepdims
    public :: hmean
    public :: hmean_keepdims
    public :: pmean
    public :: pmean_keepdims

contains

    pure subroutine reduce_assumed_rank_mean(x, axis, mean_kind, reduced, out_shape, ok, &
            p, nan_policy)
        real(dp), intent(in) :: x(..) !! observations with rank 5 to 15
        integer, intent(in) :: axis !! Fortran axis to reduce
        integer, intent(in) :: mean_kind !! internal selector for geometric, harmonic, or power mean
        real(dp), allocatable, intent(out) :: reduced(:) !! flattened reduction result
        integer, intent(out) :: out_shape(:) !! result shape with reduced axis removed
        logical, intent(out) :: ok !! true when rank and axis are valid
        real(dp), intent(in), optional :: p !! finite power for the power mean
        character(len=*), intent(in), optional :: nan_policy !! optional propagate/omit/raise policy

        integer :: base_index
        integer :: coordinate
        integer :: d
        integer, allocatable :: dims(:)
        real(dp), allocatable :: flat(:)
        integer :: i
        integer :: linear
        integer :: n_out
        integer :: n_reduce
        integer :: out_dim
        real(dp), allocatable :: slice(:)
        integer, allocatable :: strides(:)

        ok = .false.
        if (rank(x) < 5 .or. rank(x) > 15) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if
        if (size(out_shape) /= rank(x) - 1) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if
        if (axis < 1 .or. axis > rank(x)) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if
        if (mean_kind == mean_kind_pmean .and. .not. present(p)) then
            allocate(reduced(0))
            out_shape = 0
            return
        end if

        select rank (a => x)
        rank (5)
            flat = reshape(a, [size(a)])
        rank (6)
            flat = reshape(a, [size(a)])
        rank (7)
            flat = reshape(a, [size(a)])
        rank (8)
            flat = reshape(a, [size(a)])
        rank (9)
            flat = reshape(a, [size(a)])
        rank (10)
            flat = reshape(a, [size(a)])
        rank (11)
            flat = reshape(a, [size(a)])
        rank (12)
            flat = reshape(a, [size(a)])
        rank (13)
            flat = reshape(a, [size(a)])
        rank (14)
            flat = reshape(a, [size(a)])
        rank (15)
            flat = reshape(a, [size(a)])
        rank default
            allocate(reduced(0))
            out_shape = 0
            return
        end select

        dims = shape(x)
        allocate(strides(size(dims)))
        strides(1) = 1
        do d = 2, size(dims)
            strides(d) = strides(d - 1) * dims(d - 1)
        end do

        out_dim = 0
        do d = 1, size(dims)
            if (d == axis) cycle
            out_dim = out_dim + 1
            out_shape(out_dim) = dims(d)
        end do
        n_reduce = dims(axis)
        n_out = product(out_shape)
        allocate(reduced(n_out))
        allocate(slice(n_reduce))

        do i = 0, n_out - 1
            linear = i
            base_index = 1
            out_dim = 0
            do d = 1, size(dims)
                if (d == axis) cycle
                out_dim = out_dim + 1
                coordinate = mod(linear, out_shape(out_dim))
                linear = linear / out_shape(out_dim)
                base_index = base_index + coordinate * strides(d)
            end do
            do d = 1, n_reduce
                slice(d) = flat(base_index + (d - 1) * strides(axis))
            end do

            select case (mean_kind)
            case (mean_kind_gmean)
                if (present(nan_policy)) then
                    reduced(i + 1) = gmean_low( &
                        slice, nan_policy=nan_policy)
                else
                    reduced(i + 1) = gmean_low(slice)
                end if
            case (mean_kind_hmean)
                if (present(nan_policy)) then
                    reduced(i + 1) = hmean_low( &
                        slice, nan_policy=nan_policy)
                else
                    reduced(i + 1) = hmean_low(slice)
                end if
            case (mean_kind_pmean)
                if (present(nan_policy)) then
                    reduced(i + 1) = pmean_low( &
                        slice, p, nan_policy=nan_policy)
                else
                    reduced(i + 1) = pmean_low(slice, p)
                end if
            case default
                reduced = 0.0_dp
                return
            end select
        end do
        ok = .true.
    end subroutine reduce_assumed_rank_mean

    pure function gmean_rank5_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank5_assumed

    pure function gmean_rank5_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank5_assumed_nan_policy

    pure function gmean_rank5_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank5_assumed_keepdims

    pure function gmean_rank5_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank5_assumed_keepdims_nan_policy

    pure function hmean_rank5_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank5_assumed

    pure function hmean_rank5_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank5_assumed_nan_policy

    pure function hmean_rank5_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank5_assumed_keepdims

    pure function hmean_rank5_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank5_assumed_keepdims_nan_policy

    pure function pmean_rank5_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank5_assumed

    pure function pmean_rank5_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank5_assumed_nan_policy

    pure function pmean_rank5_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank5_assumed_keepdims

    pure function pmean_rank5_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank5_assumed_keepdims_nan_policy

    pure function gmean_rank6_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank6_assumed

    pure function gmean_rank6_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank6_assumed_nan_policy

    pure function gmean_rank6_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank6_assumed_keepdims

    pure function gmean_rank6_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank6_assumed_keepdims_nan_policy

    pure function hmean_rank6_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank6_assumed

    pure function hmean_rank6_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank6_assumed_nan_policy

    pure function hmean_rank6_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank6_assumed_keepdims

    pure function hmean_rank6_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank6_assumed_keepdims_nan_policy

    pure function pmean_rank6_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank6_assumed

    pure function pmean_rank6_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank6_assumed_nan_policy

    pure function pmean_rank6_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank6_assumed_keepdims

    pure function pmean_rank6_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank6_assumed_keepdims_nan_policy

    pure function gmean_rank7_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank7_assumed

    pure function gmean_rank7_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank7_assumed_nan_policy

    pure function gmean_rank7_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank7_assumed_keepdims

    pure function gmean_rank7_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank7_assumed_keepdims_nan_policy

    pure function hmean_rank7_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank7_assumed

    pure function hmean_rank7_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank7_assumed_nan_policy

    pure function hmean_rank7_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank7_assumed_keepdims

    pure function hmean_rank7_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank7_assumed_keepdims_nan_policy

    pure function pmean_rank7_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank7_assumed

    pure function pmean_rank7_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank7_assumed_nan_policy

    pure function pmean_rank7_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank7_assumed_keepdims

    pure function pmean_rank7_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank7_assumed_keepdims_nan_policy

    pure function gmean_rank8_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank8_assumed

    pure function gmean_rank8_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank8_assumed_nan_policy

    pure function gmean_rank8_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank8_assumed_keepdims

    pure function gmean_rank8_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank8_assumed_keepdims_nan_policy

    pure function hmean_rank8_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank8_assumed

    pure function hmean_rank8_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank8_assumed_nan_policy

    pure function hmean_rank8_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank8_assumed_keepdims

    pure function hmean_rank8_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank8_assumed_keepdims_nan_policy

    pure function pmean_rank8_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank8_assumed

    pure function pmean_rank8_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank8_assumed_nan_policy

    pure function pmean_rank8_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank8_assumed_keepdims

    pure function pmean_rank8_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank8_assumed_keepdims_nan_policy

    pure function gmean_rank9_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank9_assumed

    pure function gmean_rank9_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank9_assumed_nan_policy

    pure function gmean_rank9_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank9_assumed_keepdims

    pure function gmean_rank9_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank9_assumed_keepdims_nan_policy

    pure function hmean_rank9_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank9_assumed

    pure function hmean_rank9_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank9_assumed_nan_policy

    pure function hmean_rank9_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank9_assumed_keepdims

    pure function hmean_rank9_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank9_assumed_keepdims_nan_policy

    pure function pmean_rank9_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank9_assumed

    pure function pmean_rank9_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank9_assumed_nan_policy

    pure function pmean_rank9_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank9_assumed_keepdims

    pure function pmean_rank9_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(8)
        integer :: result_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank9_assumed_keepdims_nan_policy

    pure function gmean_rank10_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank10_assumed

    pure function gmean_rank10_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank10_assumed_nan_policy

    pure function gmean_rank10_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank10_assumed_keepdims

    pure function gmean_rank10_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank10_assumed_keepdims_nan_policy

    pure function hmean_rank10_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank10_assumed

    pure function hmean_rank10_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank10_assumed_nan_policy

    pure function hmean_rank10_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank10_assumed_keepdims

    pure function hmean_rank10_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank10_assumed_keepdims_nan_policy

    pure function pmean_rank10_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank10_assumed

    pure function pmean_rank10_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank10_assumed_nan_policy

    pure function pmean_rank10_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank10_assumed_keepdims

    pure function pmean_rank10_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank10_assumed_keepdims_nan_policy

    pure function gmean_rank11_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank11_assumed

    pure function gmean_rank11_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank11_assumed_nan_policy

    pure function gmean_rank11_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank11_assumed_keepdims

    pure function gmean_rank11_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank11_assumed_keepdims_nan_policy

    pure function hmean_rank11_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank11_assumed

    pure function hmean_rank11_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank11_assumed_nan_policy

    pure function hmean_rank11_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank11_assumed_keepdims

    pure function hmean_rank11_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank11_assumed_keepdims_nan_policy

    pure function pmean_rank11_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank11_assumed

    pure function pmean_rank11_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank11_assumed_nan_policy

    pure function pmean_rank11_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank11_assumed_keepdims

    pure function pmean_rank11_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank11_assumed_keepdims_nan_policy

    pure function gmean_rank12_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank12_assumed

    pure function gmean_rank12_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank12_assumed_nan_policy

    pure function gmean_rank12_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank12_assumed_keepdims

    pure function gmean_rank12_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank12_assumed_keepdims_nan_policy

    pure function hmean_rank12_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank12_assumed

    pure function hmean_rank12_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank12_assumed_nan_policy

    pure function hmean_rank12_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank12_assumed_keepdims

    pure function hmean_rank12_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank12_assumed_keepdims_nan_policy

    pure function pmean_rank12_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank12_assumed

    pure function pmean_rank12_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank12_assumed_nan_policy

    pure function pmean_rank12_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank12_assumed_keepdims

    pure function pmean_rank12_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank12_assumed_keepdims_nan_policy

    pure function gmean_rank13_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank13_assumed

    pure function gmean_rank13_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank13_assumed_nan_policy

    pure function gmean_rank13_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank13_assumed_keepdims

    pure function gmean_rank13_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank13_assumed_keepdims_nan_policy

    pure function hmean_rank13_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank13_assumed

    pure function hmean_rank13_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank13_assumed_nan_policy

    pure function hmean_rank13_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank13_assumed_keepdims

    pure function hmean_rank13_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank13_assumed_keepdims_nan_policy

    pure function pmean_rank13_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank13_assumed

    pure function pmean_rank13_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank13_assumed_nan_policy

    pure function pmean_rank13_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank13_assumed_keepdims

    pure function pmean_rank13_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank13_assumed_keepdims_nan_policy

    pure function gmean_rank14_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank14_assumed

    pure function gmean_rank14_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank14_assumed_nan_policy

    pure function gmean_rank14_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank14_assumed_keepdims

    pure function gmean_rank14_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank14_assumed_keepdims_nan_policy

    pure function hmean_rank14_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank14_assumed

    pure function hmean_rank14_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank14_assumed_nan_policy

    pure function hmean_rank14_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank14_assumed_keepdims

    pure function hmean_rank14_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank14_assumed_keepdims_nan_policy

    pure function pmean_rank14_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank14_assumed

    pure function pmean_rank14_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank14_assumed_nan_policy

    pure function pmean_rank14_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank14_assumed_keepdims

    pure function pmean_rank14_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank14_assumed_keepdims_nan_policy

    pure function gmean_rank15_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank15_assumed

    pure function gmean_rank15_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank15_assumed_nan_policy

    pure function gmean_rank15_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank15_assumed_keepdims

    pure function gmean_rank15_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_gmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank15_assumed_keepdims_nan_policy

    pure function hmean_rank15_assumed(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank15_assumed

    pure function hmean_rank15_assumed_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank15_assumed_nan_policy

    pure function hmean_rank15_assumed_keepdims(x, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank15_assumed_keepdims

    pure function hmean_rank15_assumed_keepdims_nan_policy(x, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_hmean, reduced, reduced_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank15_assumed_keepdims_nan_policy

    pure function pmean_rank15_assumed(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank15_assumed

    pure function pmean_rank15_assumed_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank15_assumed_nan_policy

    pure function pmean_rank15_assumed_keepdims(x, p, axis) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank15_assumed_keepdims

    pure function pmean_rank15_assumed_keepdims_nan_policy(x, p, axis, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! nonnegative observations; reduced along axis
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_mean(x, axis_value, mean_kind_pmean, reduced, reduced_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = shape(x)
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank15_assumed_keepdims_nan_policy

end module scifort_descriptive_highrank_means
