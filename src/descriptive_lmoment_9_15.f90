! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_lmoment_9_15
    use scifort_descriptive_lmoment_engine, only : lmoment_all_orders_assumed_rank, &
        lmoment_all_scalar_assumed_rank, lmoment_reduce_orders_assumed_rank, &
        lmoment_reduce_scalar_assumed_rank
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter :: default_orders(4) = [1, 2, 3, 4]

    interface lmoment
        module procedure lmoment_default_rank9
        module procedure lmoment_scalar_rank9
        module procedure lmoment_orders_rank9
        module procedure lmoment_default_rank10
        module procedure lmoment_scalar_rank10
        module procedure lmoment_orders_rank10
        module procedure lmoment_default_rank11
        module procedure lmoment_scalar_rank11
        module procedure lmoment_orders_rank11
        module procedure lmoment_default_rank12
        module procedure lmoment_scalar_rank12
        module procedure lmoment_orders_rank12
        module procedure lmoment_default_rank13
        module procedure lmoment_scalar_rank13
        module procedure lmoment_orders_rank13
        module procedure lmoment_default_rank14
        module procedure lmoment_scalar_rank14
        module procedure lmoment_orders_rank14
        module procedure lmoment_default_rank15
        module procedure lmoment_scalar_rank15
        module procedure lmoment_orders_rank15
    end interface lmoment

    interface lmoment_all
        module procedure lmoment_all_default_rank9
        module procedure lmoment_all_scalar_rank9
        module procedure lmoment_all_orders_rank9
        module procedure lmoment_all_default_rank10
        module procedure lmoment_all_scalar_rank10
        module procedure lmoment_all_orders_rank10
        module procedure lmoment_all_default_rank11
        module procedure lmoment_all_scalar_rank11
        module procedure lmoment_all_orders_rank11
        module procedure lmoment_all_default_rank12
        module procedure lmoment_all_scalar_rank12
        module procedure lmoment_all_orders_rank12
        module procedure lmoment_all_default_rank13
        module procedure lmoment_all_scalar_rank13
        module procedure lmoment_all_orders_rank13
        module procedure lmoment_all_default_rank14
        module procedure lmoment_all_scalar_rank14
        module procedure lmoment_all_orders_rank14
        module procedure lmoment_all_default_rank15
        module procedure lmoment_all_scalar_rank15
        module procedure lmoment_all_orders_rank15
    end interface lmoment_all

    interface lmoment_keepdims
        module procedure lmoment_keepdims_scalar_rank9
        module procedure lmoment_keepdims_default_rank9
        module procedure lmoment_keepdims_orders_rank9
        module procedure lmoment_keepdims_scalar_rank10
        module procedure lmoment_keepdims_default_rank10
        module procedure lmoment_keepdims_orders_rank10
        module procedure lmoment_keepdims_scalar_rank11
        module procedure lmoment_keepdims_default_rank11
        module procedure lmoment_keepdims_orders_rank11
        module procedure lmoment_keepdims_scalar_rank12
        module procedure lmoment_keepdims_default_rank12
        module procedure lmoment_keepdims_orders_rank12
        module procedure lmoment_keepdims_scalar_rank13
        module procedure lmoment_keepdims_default_rank13
        module procedure lmoment_keepdims_orders_rank13
        module procedure lmoment_keepdims_scalar_rank14
        module procedure lmoment_keepdims_default_rank14
        module procedure lmoment_keepdims_orders_rank14
        module procedure lmoment_keepdims_scalar_rank15
    end interface lmoment_keepdims

    public :: lmoment
    public :: lmoment_all
    public :: lmoment_keepdims

contains

    pure function lmoment_default_rank9(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(9)
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:9) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank9

    pure function lmoment_scalar_rank9(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank9

    pure function lmoment_orders_rank9(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(9)
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:9) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank9

    pure function lmoment_all_default_rank9(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank9

    pure function lmoment_all_scalar_rank9(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank9

    pure function lmoment_all_orders_rank9(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank9

    pure function lmoment_keepdims_default_rank9(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(9)
        integer :: full_shape(10)
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:9) = out_shape
        j = 1
        do d = 1, 10
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank9

    pure function lmoment_keepdims_scalar_rank9(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(9)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 9
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank9

    pure function lmoment_keepdims_orders_rank9(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(10)
        integer :: j
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 9
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank9

    pure function lmoment_default_rank10(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(10)
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:10) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank10

    pure function lmoment_scalar_rank10(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank10

    pure function lmoment_orders_rank10(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(10)
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:10) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank10

    pure function lmoment_all_default_rank10(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank10

    pure function lmoment_all_scalar_rank10(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank10

    pure function lmoment_all_orders_rank10(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank10

    pure function lmoment_keepdims_default_rank10(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(10)
        integer :: full_shape(11)
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:10) = out_shape
        j = 1
        do d = 1, 11
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank10

    pure function lmoment_keepdims_scalar_rank10(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(10)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 10
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank10

    pure function lmoment_keepdims_orders_rank10(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(11)
        integer :: j
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 10
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank10

    pure function lmoment_default_rank11(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(11)
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:11) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank11

    pure function lmoment_scalar_rank11(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank11

    pure function lmoment_orders_rank11(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(11)
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:11) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank11

    pure function lmoment_all_default_rank11(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank11

    pure function lmoment_all_scalar_rank11(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank11

    pure function lmoment_all_orders_rank11(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank11

    pure function lmoment_keepdims_default_rank11(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(11)
        integer :: full_shape(12)
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:11) = out_shape
        j = 1
        do d = 1, 12
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank11

    pure function lmoment_keepdims_scalar_rank11(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(11)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 11
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank11

    pure function lmoment_keepdims_orders_rank11(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(12)
        integer :: j
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 11
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank11

    pure function lmoment_default_rank12(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(12)
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:12) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank12

    pure function lmoment_scalar_rank12(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank12

    pure function lmoment_orders_rank12(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(12)
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:12) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank12

    pure function lmoment_all_default_rank12(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank12

    pure function lmoment_all_scalar_rank12(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank12

    pure function lmoment_all_orders_rank12(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank12

    pure function lmoment_keepdims_default_rank12(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(12)
        integer :: full_shape(13)
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:12) = out_shape
        j = 1
        do d = 1, 13
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank12

    pure function lmoment_keepdims_scalar_rank12(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(12)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 12
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank12

    pure function lmoment_keepdims_orders_rank12(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(13)
        integer :: j
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 12
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank12

    pure function lmoment_default_rank13(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(13)
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:13) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank13

    pure function lmoment_scalar_rank13(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank13

    pure function lmoment_orders_rank13(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(13)
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:13) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank13

    pure function lmoment_all_default_rank13(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank13

    pure function lmoment_all_scalar_rank13(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank13

    pure function lmoment_all_orders_rank13(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank13

    pure function lmoment_keepdims_default_rank13(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(13)
        integer :: full_shape(14)
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:13) = out_shape
        j = 1
        do d = 1, 14
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank13

    pure function lmoment_keepdims_scalar_rank13(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(13)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 13
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank13

    pure function lmoment_keepdims_orders_rank13(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(14)
        integer :: j
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 13
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank13

    pure function lmoment_default_rank14(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(14)
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:14) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank14

    pure function lmoment_scalar_rank14(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank14

    pure function lmoment_orders_rank14(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(14)
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:14) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank14

    pure function lmoment_all_default_rank14(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank14

    pure function lmoment_all_scalar_rank14(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank14

    pure function lmoment_all_orders_rank14(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank14

    pure function lmoment_keepdims_default_rank14(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(14)
        integer :: full_shape(15)
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:14) = out_shape
        j = 1
        do d = 1, 15
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank14

    pure function lmoment_keepdims_scalar_rank14(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(14)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 14
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank14

    pure function lmoment_keepdims_orders_rank14(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(15)
        integer :: j
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 14
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank14

    pure function lmoment_default_rank15(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(15)
        logical :: ok
        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:15) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank15

    pure function lmoment_scalar_rank15(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank15

    pure function lmoment_orders_rank15(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(15)
        logical :: ok
        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:15) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank15

    pure function lmoment_all_default_rank15(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank15

    pure function lmoment_all_scalar_rank15(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank15

    pure function lmoment_all_orders_rank15(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank15

    pure function lmoment_keepdims_scalar_rank15(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(15)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 15
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank15

end module scifort_descriptive_lmoment_9_15
