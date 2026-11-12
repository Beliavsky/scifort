! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_highrank_weights_4_9
    use scifort_descriptive_weighted_engine, only : mean_kind_gmean, mean_kind_hmean, &
        mean_kind_pmean, reduce_assumed_rank_weighted_mean
    use scifort_kinds, only : dp
    implicit none
    private

    interface gmean
        module procedure gmean_rank4_weighted
        module procedure gmean_rank4_weighted_nan_policy
        module procedure gmean_rank5_weighted
        module procedure gmean_rank5_weighted_nan_policy
        module procedure gmean_rank6_weighted
        module procedure gmean_rank6_weighted_nan_policy
        module procedure gmean_rank7_weighted
        module procedure gmean_rank7_weighted_nan_policy
        module procedure gmean_rank8_weighted
        module procedure gmean_rank8_weighted_nan_policy
        module procedure gmean_rank9_weighted
        module procedure gmean_rank9_weighted_nan_policy
    end interface gmean

    interface hmean
        module procedure hmean_rank4_weighted
        module procedure hmean_rank4_weighted_nan_policy
        module procedure hmean_rank5_weighted
        module procedure hmean_rank5_weighted_nan_policy
        module procedure hmean_rank6_weighted
        module procedure hmean_rank6_weighted_nan_policy
        module procedure hmean_rank7_weighted
        module procedure hmean_rank7_weighted_nan_policy
        module procedure hmean_rank8_weighted
        module procedure hmean_rank8_weighted_nan_policy
        module procedure hmean_rank9_weighted
        module procedure hmean_rank9_weighted_nan_policy
    end interface hmean

    interface pmean
        module procedure pmean_rank4_weighted
        module procedure pmean_rank4_weighted_nan_policy
        module procedure pmean_rank5_weighted
        module procedure pmean_rank5_weighted_nan_policy
        module procedure pmean_rank6_weighted
        module procedure pmean_rank6_weighted_nan_policy
        module procedure pmean_rank7_weighted
        module procedure pmean_rank7_weighted_nan_policy
        module procedure pmean_rank8_weighted
        module procedure pmean_rank8_weighted_nan_policy
        module procedure pmean_rank9_weighted
        module procedure pmean_rank9_weighted_nan_policy
    end interface pmean

    interface gmean_keepdims
        module procedure gmean_rank4_weighted_keepdims
        module procedure gmean_rank4_weighted_keepdims_nan_policy
        module procedure gmean_rank5_weighted_keepdims
        module procedure gmean_rank5_weighted_keepdims_nan_policy
        module procedure gmean_rank6_weighted_keepdims
        module procedure gmean_rank6_weighted_keepdims_nan_policy
        module procedure gmean_rank7_weighted_keepdims
        module procedure gmean_rank7_weighted_keepdims_nan_policy
        module procedure gmean_rank8_weighted_keepdims
        module procedure gmean_rank8_weighted_keepdims_nan_policy
        module procedure gmean_rank9_weighted_keepdims
        module procedure gmean_rank9_weighted_keepdims_nan_policy
    end interface gmean_keepdims

    interface hmean_keepdims
        module procedure hmean_rank4_weighted_keepdims
        module procedure hmean_rank4_weighted_keepdims_nan_policy
        module procedure hmean_rank5_weighted_keepdims
        module procedure hmean_rank5_weighted_keepdims_nan_policy
        module procedure hmean_rank6_weighted_keepdims
        module procedure hmean_rank6_weighted_keepdims_nan_policy
        module procedure hmean_rank7_weighted_keepdims
        module procedure hmean_rank7_weighted_keepdims_nan_policy
        module procedure hmean_rank8_weighted_keepdims
        module procedure hmean_rank8_weighted_keepdims_nan_policy
        module procedure hmean_rank9_weighted_keepdims
        module procedure hmean_rank9_weighted_keepdims_nan_policy
    end interface hmean_keepdims

    interface pmean_keepdims
        module procedure pmean_rank4_weighted_keepdims
        module procedure pmean_rank4_weighted_keepdims_nan_policy
        module procedure pmean_rank5_weighted_keepdims
        module procedure pmean_rank5_weighted_keepdims_nan_policy
        module procedure pmean_rank6_weighted_keepdims
        module procedure pmean_rank6_weighted_keepdims_nan_policy
        module procedure pmean_rank7_weighted_keepdims
        module procedure pmean_rank7_weighted_keepdims_nan_policy
        module procedure pmean_rank8_weighted_keepdims
        module procedure pmean_rank8_weighted_keepdims_nan_policy
        module procedure pmean_rank9_weighted_keepdims
        module procedure pmean_rank9_weighted_keepdims_nan_policy
    end interface pmean_keepdims

    public :: gmean
    public :: hmean
    public :: pmean
    public :: gmean_keepdims
    public :: hmean_keepdims
    public :: pmean_keepdims

contains

    pure function gmean_rank4_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank4_weighted

    pure function gmean_rank4_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank4_weighted_nan_policy

    pure function gmean_rank4_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        integer :: result_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank4_weighted_keepdims

    pure function gmean_rank4_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        integer :: result_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank4_weighted_keepdims_nan_policy

    pure function hmean_rank4_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank4_weighted

    pure function hmean_rank4_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank4_weighted_nan_policy

    pure function hmean_rank4_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        integer :: result_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank4_weighted_keepdims

    pure function hmean_rank4_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        integer :: result_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank4_weighted_keepdims_nan_policy

    pure function pmean_rank4_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank4_weighted

    pure function pmean_rank4_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank4_weighted_nan_policy

    pure function pmean_rank4_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        integer :: result_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank4_weighted_keepdims

    pure function pmean_rank4_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 4; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(4)
        integer :: reduced_shape(3)
        integer :: result_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank4_weighted_keepdims_nan_policy

    pure function gmean_rank5_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank5_weighted

    pure function gmean_rank5_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank5_weighted_nan_policy

    pure function gmean_rank5_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank5_weighted_keepdims

    pure function gmean_rank5_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank5_weighted_keepdims_nan_policy

    pure function hmean_rank5_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank5_weighted

    pure function hmean_rank5_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank5_weighted_nan_policy

    pure function hmean_rank5_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank5_weighted_keepdims

    pure function hmean_rank5_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank5_weighted_keepdims_nan_policy

    pure function pmean_rank5_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank5_weighted

    pure function pmean_rank5_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank5_weighted_nan_policy

    pure function pmean_rank5_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank5_weighted_keepdims

    pure function pmean_rank5_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 5; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(5)
        integer :: reduced_shape(4)
        integer :: result_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank5_weighted_keepdims_nan_policy

    pure function gmean_rank6_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank6_weighted

    pure function gmean_rank6_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank6_weighted_nan_policy

    pure function gmean_rank6_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank6_weighted_keepdims

    pure function gmean_rank6_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank6_weighted_keepdims_nan_policy

    pure function hmean_rank6_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank6_weighted

    pure function hmean_rank6_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank6_weighted_nan_policy

    pure function hmean_rank6_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank6_weighted_keepdims

    pure function hmean_rank6_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank6_weighted_keepdims_nan_policy

    pure function pmean_rank6_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank6_weighted

    pure function pmean_rank6_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank6_weighted_nan_policy

    pure function pmean_rank6_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank6_weighted_keepdims

    pure function pmean_rank6_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 6; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(6)
        integer :: reduced_shape(5)
        integer :: result_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank6_weighted_keepdims_nan_policy

    pure function gmean_rank7_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank7_weighted

    pure function gmean_rank7_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank7_weighted_nan_policy

    pure function gmean_rank7_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank7_weighted_keepdims

    pure function gmean_rank7_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank7_weighted_keepdims_nan_policy

    pure function hmean_rank7_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank7_weighted

    pure function hmean_rank7_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank7_weighted_nan_policy

    pure function hmean_rank7_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank7_weighted_keepdims

    pure function hmean_rank7_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank7_weighted_keepdims_nan_policy

    pure function pmean_rank7_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank7_weighted

    pure function pmean_rank7_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank7_weighted_nan_policy

    pure function pmean_rank7_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank7_weighted_keepdims

    pure function pmean_rank7_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 7; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(7)
        integer :: reduced_shape(6)
        integer :: result_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank7_weighted_keepdims_nan_policy

    pure function gmean_rank8_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank8_weighted

    pure function gmean_rank8_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank8_weighted_nan_policy

    pure function gmean_rank8_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank8_weighted_keepdims

    pure function gmean_rank8_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank8_weighted_keepdims_nan_policy

    pure function hmean_rank8_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank8_weighted

    pure function hmean_rank8_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank8_weighted_nan_policy

    pure function hmean_rank8_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank8_weighted_keepdims

    pure function hmean_rank8_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank8_weighted_keepdims_nan_policy

    pure function pmean_rank8_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank8_weighted

    pure function pmean_rank8_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank8_weighted_nan_policy

    pure function pmean_rank8_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank8_weighted_keepdims

    pure function pmean_rank8_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 8; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(8)
        integer :: reduced_shape(7)
        integer :: result_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank8_weighted_keepdims_nan_policy

    pure function gmean_rank9_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank9_weighted

    pure function gmean_rank9_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_gmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function gmean_rank9_weighted_nan_policy

    pure function gmean_rank9_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        integer :: result_shape(9)
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
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank9_weighted_keepdims

    pure function gmean_rank9_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        integer :: result_shape(9)
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
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function gmean_rank9_weighted_keepdims_nan_policy

    pure function hmean_rank9_weighted(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank9_weighted

    pure function hmean_rank9_weighted_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_hmean, reduced, reduced_shape, broadcast_shape, ok, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function hmean_rank9_weighted_nan_policy

    pure function hmean_rank9_weighted_keepdims(x, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        integer :: result_shape(9)
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
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank9_weighted_keepdims

    pure function hmean_rank9_weighted_keepdims_nan_policy(x, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        integer :: result_shape(9)
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
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function hmean_rank9_weighted_keepdims_nan_policy

    pure function pmean_rank9_weighted(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank9_weighted

    pure function pmean_rank9_weighted_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        logical :: ok
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call reduce_assumed_rank_weighted_mean(x, weights, axis_value, &
            mean_kind_pmean, reduced, reduced_shape, broadcast_shape, ok, p=p, nan_policy=nan_policy)
        if (.not. ok) then
            allocate(values(0, 0, 0, 0, 0, 0, 0, 0))
            return
        end if
        values = reshape(reduced, reduced_shape)
    end function pmean_rank9_weighted_nan_policy

    pure function pmean_rank9_weighted_keepdims(x, p, axis, weights) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        integer :: result_shape(9)
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
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank9_weighted_keepdims

    pure function pmean_rank9_weighted_keepdims_nan_policy(x, p, axis, weights, nan_policy) result(values)
        real(dp), intent(in) :: x(:, :, :, :, :, :, :, :, :) !! observations
        real(dp), intent(in) :: p !! finite power
        integer, intent(in), optional :: axis !! Fortran axis, 1 to 9; default 1
        real(dp), intent(in) :: weights(..) !! scalar or trailing-broadcast weights
        character(len=*), intent(in) :: nan_policy !! 'propagate', 'omit', or 'raise'
        real(dp), allocatable :: values(:, :, :, :, :, :, :, :, :)

        integer :: axis_value
        integer :: broadcast_shape(9)
        integer :: reduced_shape(8)
        integer :: result_shape(9)
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
        result_shape = broadcast_shape
        result_shape(axis_value) = 1
        values = reshape(reduced, result_shape)
    end function pmean_rank9_weighted_keepdims_nan_policy

end module scifort_descriptive_highrank_weights_4_9
