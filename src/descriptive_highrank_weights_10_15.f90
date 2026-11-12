! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_highrank_weights_10_15
    use scifort_descriptive_weighted_engine, only : mean_kind_gmean, mean_kind_hmean, &
        mean_kind_pmean, reduce_assumed_rank_weighted_mean
    use scifort_kinds, only : dp
    implicit none
    private

    interface gmean
        module procedure gmean_rank10_weighted
        module procedure gmean_rank10_weighted_nan_policy
        module procedure gmean_rank11_weighted
        module procedure gmean_rank11_weighted_nan_policy
        module procedure gmean_rank12_weighted
        module procedure gmean_rank12_weighted_nan_policy
        module procedure gmean_rank13_weighted
        module procedure gmean_rank13_weighted_nan_policy
        module procedure gmean_rank14_weighted
        module procedure gmean_rank14_weighted_nan_policy
        module procedure gmean_rank15_weighted
        module procedure gmean_rank15_weighted_nan_policy
    end interface gmean

    interface hmean
        module procedure hmean_rank10_weighted
        module procedure hmean_rank10_weighted_nan_policy
        module procedure hmean_rank11_weighted
        module procedure hmean_rank11_weighted_nan_policy
        module procedure hmean_rank12_weighted
        module procedure hmean_rank12_weighted_nan_policy
        module procedure hmean_rank13_weighted
        module procedure hmean_rank13_weighted_nan_policy
        module procedure hmean_rank14_weighted
        module procedure hmean_rank14_weighted_nan_policy
        module procedure hmean_rank15_weighted
        module procedure hmean_rank15_weighted_nan_policy
    end interface hmean

    interface pmean
        module procedure pmean_rank10_weighted
        module procedure pmean_rank10_weighted_nan_policy
        module procedure pmean_rank11_weighted
        module procedure pmean_rank11_weighted_nan_policy
        module procedure pmean_rank12_weighted
        module procedure pmean_rank12_weighted_nan_policy
        module procedure pmean_rank13_weighted
        module procedure pmean_rank13_weighted_nan_policy
        module procedure pmean_rank14_weighted
        module procedure pmean_rank14_weighted_nan_policy
        module procedure pmean_rank15_weighted
        module procedure pmean_rank15_weighted_nan_policy
    end interface pmean

    interface gmean_keepdims
        module procedure gmean_rank10_weighted_keepdims
        module procedure gmean_rank10_weighted_keepdims_nan_policy
        module procedure gmean_rank11_weighted_keepdims
        module procedure gmean_rank11_weighted_keepdims_nan_policy
        module procedure gmean_rank12_weighted_keepdims
        module procedure gmean_rank12_weighted_keepdims_nan_policy
        module procedure gmean_rank13_weighted_keepdims
        module procedure gmean_rank13_weighted_keepdims_nan_policy
        module procedure gmean_rank14_weighted_keepdims
        module procedure gmean_rank14_weighted_keepdims_nan_policy
        module procedure gmean_rank15_weighted_keepdims
        module procedure gmean_rank15_weighted_keepdims_nan_policy
    end interface gmean_keepdims

    interface hmean_keepdims
        module procedure hmean_rank10_weighted_keepdims
        module procedure hmean_rank10_weighted_keepdims_nan_policy
        module procedure hmean_rank11_weighted_keepdims
        module procedure hmean_rank11_weighted_keepdims_nan_policy
        module procedure hmean_rank12_weighted_keepdims
        module procedure hmean_rank12_weighted_keepdims_nan_policy
        module procedure hmean_rank13_weighted_keepdims
        module procedure hmean_rank13_weighted_keepdims_nan_policy
        module procedure hmean_rank14_weighted_keepdims
        module procedure hmean_rank14_weighted_keepdims_nan_policy
        module procedure hmean_rank15_weighted_keepdims
        module procedure hmean_rank15_weighted_keepdims_nan_policy
    end interface hmean_keepdims

    interface pmean_keepdims
        module procedure pmean_rank10_weighted_keepdims
        module procedure pmean_rank10_weighted_keepdims_nan_policy
        module procedure pmean_rank11_weighted_keepdims
        module procedure pmean_rank11_weighted_keepdims_nan_policy
        module procedure pmean_rank12_weighted_keepdims
        module procedure pmean_rank12_weighted_keepdims_nan_policy
        module procedure pmean_rank13_weighted_keepdims
        module procedure pmean_rank13_weighted_keepdims_nan_policy
        module procedure pmean_rank14_weighted_keepdims
        module procedure pmean_rank14_weighted_keepdims_nan_policy
        module procedure pmean_rank15_weighted_keepdims
        module procedure pmean_rank15_weighted_keepdims_nan_policy
    end interface pmean_keepdims

    public :: gmean
    public :: hmean
    public :: pmean
    public :: gmean_keepdims
    public :: hmean_keepdims
    public :: pmean_keepdims

contains

    pure function gmean_rank10_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank10_weighted

    pure function gmean_rank10_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank10_weighted_nan_policy

    pure function gmean_rank10_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank10_weighted_keepdims

    pure function gmean_rank10_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank10_weighted_keepdims_nan_policy

    pure function hmean_rank10_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank10_weighted

    pure function hmean_rank10_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank10_weighted_nan_policy

    pure function hmean_rank10_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank10_weighted_keepdims

    pure function hmean_rank10_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank10_weighted_keepdims_nan_policy

    pure function pmean_rank10_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank10_weighted

    pure function pmean_rank10_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank10_weighted_nan_policy

    pure function pmean_rank10_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank10_weighted_keepdims

    pure function pmean_rank10_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 10; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(10)
        integer :: reduced_shape(9)
        integer :: result_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank10_weighted_keepdims_nan_policy

    pure function gmean_rank11_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank11_weighted

    pure function gmean_rank11_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank11_weighted_nan_policy

    pure function gmean_rank11_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank11_weighted_keepdims

    pure function gmean_rank11_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank11_weighted_keepdims_nan_policy

    pure function hmean_rank11_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank11_weighted

    pure function hmean_rank11_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank11_weighted_nan_policy

    pure function hmean_rank11_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank11_weighted_keepdims

    pure function hmean_rank11_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank11_weighted_keepdims_nan_policy

    pure function pmean_rank11_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank11_weighted

    pure function pmean_rank11_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank11_weighted_nan_policy

    pure function pmean_rank11_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank11_weighted_keepdims

    pure function pmean_rank11_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 11; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(11)
        integer :: reduced_shape(10)
        integer :: result_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank11_weighted_keepdims_nan_policy

    pure function gmean_rank12_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank12_weighted

    pure function gmean_rank12_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank12_weighted_nan_policy

    pure function gmean_rank12_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank12_weighted_keepdims

    pure function gmean_rank12_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank12_weighted_keepdims_nan_policy

    pure function hmean_rank12_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank12_weighted

    pure function hmean_rank12_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank12_weighted_nan_policy

    pure function hmean_rank12_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank12_weighted_keepdims

    pure function hmean_rank12_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank12_weighted_keepdims_nan_policy

    pure function pmean_rank12_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank12_weighted

    pure function pmean_rank12_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank12_weighted_nan_policy

    pure function pmean_rank12_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank12_weighted_keepdims

    pure function pmean_rank12_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 12; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(12)
        integer :: reduced_shape(11)
        integer :: result_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank12_weighted_keepdims_nan_policy

    pure function gmean_rank13_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank13_weighted

    pure function gmean_rank13_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank13_weighted_nan_policy

    pure function gmean_rank13_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank13_weighted_keepdims

    pure function gmean_rank13_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank13_weighted_keepdims_nan_policy

    pure function hmean_rank13_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank13_weighted

    pure function hmean_rank13_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank13_weighted_nan_policy

    pure function hmean_rank13_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank13_weighted_keepdims

    pure function hmean_rank13_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank13_weighted_keepdims_nan_policy

    pure function pmean_rank13_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank13_weighted

    pure function pmean_rank13_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank13_weighted_nan_policy

    pure function pmean_rank13_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank13_weighted_keepdims

    pure function pmean_rank13_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 13; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(13)
        integer :: reduced_shape(12)
        integer :: result_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank13_weighted_keepdims_nan_policy

    pure function gmean_rank14_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank14_weighted

    pure function gmean_rank14_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank14_weighted_nan_policy

    pure function gmean_rank14_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank14_weighted_keepdims

    pure function gmean_rank14_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank14_weighted_keepdims_nan_policy

    pure function hmean_rank14_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank14_weighted

    pure function hmean_rank14_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank14_weighted_nan_policy

    pure function hmean_rank14_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank14_weighted_keepdims

    pure function hmean_rank14_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank14_weighted_keepdims_nan_policy

    pure function pmean_rank14_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank14_weighted

    pure function pmean_rank14_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank14_weighted_nan_policy

    pure function pmean_rank14_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank14_weighted_keepdims

    pure function pmean_rank14_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 14; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(14)
        integer :: reduced_shape(13)
        integer :: result_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank14_weighted_keepdims_nan_policy

    pure function gmean_rank15_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank15_weighted

    pure function gmean_rank15_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank15_weighted_nan_policy

    pure function gmean_rank15_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank15_weighted_keepdims

    pure function gmean_rank15_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank15_weighted_keepdims_nan_policy

    pure function hmean_rank15_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank15_weighted

    pure function hmean_rank15_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank15_weighted_nan_policy

    pure function hmean_rank15_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank15_weighted_keepdims

    pure function hmean_rank15_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank15_weighted_keepdims_nan_policy

    pure function pmean_rank15_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank15_weighted

    pure function pmean_rank15_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank15_weighted_nan_policy

    pure function pmean_rank15_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank15_weighted_keepdims

    pure function pmean_rank15_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 15; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(15)
        integer :: reduced_shape(14)
        integer :: result_shape(15)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank15_weighted_keepdims_nan_policy

end module scifort_descriptive_highrank_weights_10_15
