! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_lmoment_1_8
    use scifort_descriptive_lmoment_engine, only : lmoment_all_orders_assumed_rank, &
        lmoment_all_scalar_assumed_rank, lmoment_reduce_orders_assumed_rank, &
        lmoment_reduce_scalar_assumed_rank
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    integer, parameter :: default_orders(4) = [1, 2, 3, 4]

    interface lmoment
        module procedure lmoment_default_rank1
        module procedure lmoment_scalar_rank1
        module procedure lmoment_orders_rank1
        module procedure lmoment_default_rank2
        module procedure lmoment_scalar_rank2
        module procedure lmoment_orders_rank2
        module procedure lmoment_default_rank3
        module procedure lmoment_scalar_rank3
        module procedure lmoment_orders_rank3
        module procedure lmoment_default_rank4
        module procedure lmoment_scalar_rank4
        module procedure lmoment_orders_rank4
        module procedure lmoment_default_rank5
        module procedure lmoment_scalar_rank5
        module procedure lmoment_orders_rank5
        module procedure lmoment_default_rank6
        module procedure lmoment_scalar_rank6
        module procedure lmoment_orders_rank6
        module procedure lmoment_default_rank7
        module procedure lmoment_scalar_rank7
        module procedure lmoment_orders_rank7
        module procedure lmoment_default_rank8
        module procedure lmoment_scalar_rank8
        module procedure lmoment_orders_rank8
    end interface lmoment

    interface lmoment_all
        module procedure lmoment_all_default_rank1
        module procedure lmoment_all_scalar_rank1
        module procedure lmoment_all_orders_rank1
        module procedure lmoment_all_default_rank2
        module procedure lmoment_all_scalar_rank2
        module procedure lmoment_all_orders_rank2
        module procedure lmoment_all_default_rank3
        module procedure lmoment_all_scalar_rank3
        module procedure lmoment_all_orders_rank3
        module procedure lmoment_all_default_rank4
        module procedure lmoment_all_scalar_rank4
        module procedure lmoment_all_orders_rank4
        module procedure lmoment_all_default_rank5
        module procedure lmoment_all_scalar_rank5
        module procedure lmoment_all_orders_rank5
        module procedure lmoment_all_default_rank6
        module procedure lmoment_all_scalar_rank6
        module procedure lmoment_all_orders_rank6
        module procedure lmoment_all_default_rank7
        module procedure lmoment_all_scalar_rank7
        module procedure lmoment_all_orders_rank7
        module procedure lmoment_all_default_rank8
        module procedure lmoment_all_scalar_rank8
        module procedure lmoment_all_orders_rank8
    end interface lmoment_all

    interface lmoment_keepdims
        module procedure lmoment_keepdims_scalar_rank1
        module procedure lmoment_keepdims_default_rank1
        module procedure lmoment_keepdims_orders_rank1
        module procedure lmoment_keepdims_scalar_rank2
        module procedure lmoment_keepdims_default_rank2
        module procedure lmoment_keepdims_orders_rank2
        module procedure lmoment_keepdims_scalar_rank3
        module procedure lmoment_keepdims_default_rank3
        module procedure lmoment_keepdims_orders_rank3
        module procedure lmoment_keepdims_scalar_rank4
        module procedure lmoment_keepdims_default_rank4
        module procedure lmoment_keepdims_orders_rank4
        module procedure lmoment_keepdims_scalar_rank5
        module procedure lmoment_keepdims_default_rank5
        module procedure lmoment_keepdims_orders_rank5
        module procedure lmoment_keepdims_scalar_rank6
        module procedure lmoment_keepdims_default_rank6
        module procedure lmoment_keepdims_orders_rank6
        module procedure lmoment_keepdims_scalar_rank7
        module procedure lmoment_keepdims_default_rank7
        module procedure lmoment_keepdims_orders_rank7
        module procedure lmoment_keepdims_scalar_rank8
        module procedure lmoment_keepdims_default_rank8
        module procedure lmoment_keepdims_orders_rank8
    end interface lmoment_keepdims

    public :: lmoment
    public :: lmoment_all
    public :: lmoment_keepdims

contains

    pure function lmoment_default_rank1(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        full_shape(1) = size(default_orders)
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank1

    pure function lmoment_scalar_rank1(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        integer :: axis_value
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            values = quiet_nan(0.0_dp)
            return
        end if
        values = reduced(1)
    end function lmoment_scalar_rank1

    pure function lmoment_orders_rank1(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        full_shape(1) = size(order)
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank1

    pure function lmoment_all_default_rank1(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank1

    pure function lmoment_all_scalar_rank1(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank1

    pure function lmoment_all_orders_rank1(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank1

    pure function lmoment_keepdims_default_rank1(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        integer :: base_shape(1)
        integer :: full_shape(2)
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        j = 1
        do d = 1, 2
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank1

    pure function lmoment_keepdims_scalar_rank1(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: full_shape(1)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        j = 1
        do d = 1, 1
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank1

    pure function lmoment_keepdims_orders_rank1(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(2)
        integer :: j
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 1
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank1

    pure function lmoment_default_rank2(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        integer :: full_shape(2)
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:2) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank2

    pure function lmoment_scalar_rank2(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank2

    pure function lmoment_orders_rank2(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        integer :: full_shape(2)
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:2) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank2

    pure function lmoment_all_default_rank2(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank2

    pure function lmoment_all_scalar_rank2(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank2

    pure function lmoment_all_orders_rank2(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank2

    pure function lmoment_keepdims_default_rank2(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        integer :: base_shape(2)
        integer :: full_shape(3)
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:2) = out_shape
        j = 1
        do d = 1, 3
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank2

    pure function lmoment_keepdims_scalar_rank2(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        integer :: full_shape(2)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        j = 1
        do d = 1, 2
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank2

    pure function lmoment_keepdims_orders_rank2(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(3)
        integer :: j
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 2
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank2

    pure function lmoment_default_rank3(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        integer :: full_shape(3)
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:3) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank3

    pure function lmoment_scalar_rank3(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank3

    pure function lmoment_orders_rank3(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        integer :: full_shape(3)
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:3) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank3

    pure function lmoment_all_default_rank3(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank3

    pure function lmoment_all_scalar_rank3(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank3

    pure function lmoment_all_orders_rank3(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank3

    pure function lmoment_keepdims_default_rank3(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        integer :: base_shape(3)
        integer :: full_shape(4)
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:3) = out_shape
        j = 1
        do d = 1, 4
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank3

    pure function lmoment_keepdims_scalar_rank3(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        integer :: full_shape(3)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        j = 1
        do d = 1, 3
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank3

    pure function lmoment_keepdims_orders_rank3(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(4)
        integer :: j
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 3
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank3

    pure function lmoment_default_rank4(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        integer :: full_shape(4)
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:4) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank4

    pure function lmoment_scalar_rank4(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank4

    pure function lmoment_orders_rank4(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        integer :: full_shape(4)
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:4) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank4

    pure function lmoment_all_default_rank4(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank4

    pure function lmoment_all_scalar_rank4(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank4

    pure function lmoment_all_orders_rank4(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank4

    pure function lmoment_keepdims_default_rank4(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(4)
        integer :: full_shape(5)
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:4) = out_shape
        j = 1
        do d = 1, 5
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank4

    pure function lmoment_keepdims_scalar_rank4(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        integer :: full_shape(4)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 4
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank4

    pure function lmoment_keepdims_orders_rank4(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(5)
        integer :: j
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 4
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank4

    pure function lmoment_default_rank5(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(5)
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:5) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank5

    pure function lmoment_scalar_rank5(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank5

    pure function lmoment_orders_rank5(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(5)
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:5) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank5

    pure function lmoment_all_default_rank5(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank5

    pure function lmoment_all_scalar_rank5(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank5

    pure function lmoment_all_orders_rank5(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank5

    pure function lmoment_keepdims_default_rank5(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(5)
        integer :: full_shape(6)
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:5) = out_shape
        j = 1
        do d = 1, 6
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank5

    pure function lmoment_keepdims_scalar_rank5(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(5)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 5
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank5

    pure function lmoment_keepdims_orders_rank5(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(6)
        integer :: j
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 5
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank5

    pure function lmoment_default_rank6(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(6)
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:6) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank6

    pure function lmoment_scalar_rank6(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank6

    pure function lmoment_orders_rank6(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(6)
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:6) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank6

    pure function lmoment_all_default_rank6(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank6

    pure function lmoment_all_scalar_rank6(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank6

    pure function lmoment_all_orders_rank6(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank6

    pure function lmoment_keepdims_default_rank6(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(6)
        integer :: full_shape(7)
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:6) = out_shape
        j = 1
        do d = 1, 7
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank6

    pure function lmoment_keepdims_scalar_rank6(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(6)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 6
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank6

    pure function lmoment_keepdims_orders_rank6(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(7)
        integer :: j
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 6
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank6

    pure function lmoment_default_rank7(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(7)
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:7) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank7

    pure function lmoment_scalar_rank7(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank7

    pure function lmoment_orders_rank7(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(7)
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:7) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank7

    pure function lmoment_all_default_rank7(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank7

    pure function lmoment_all_scalar_rank7(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank7

    pure function lmoment_all_orders_rank7(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank7

    pure function lmoment_keepdims_default_rank7(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(7)
        integer :: full_shape(8)
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:7) = out_shape
        j = 1
        do d = 1, 8
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank7

    pure function lmoment_keepdims_scalar_rank7(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(7)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 7
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank7

    pure function lmoment_keepdims_orders_rank7(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(8)
        integer :: j
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        j = 1
        do d = 1, 7
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank7

    pure function lmoment_default_rank8(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(8)
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(default_orders)
        full_shape(2:8) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_default_rank8

    pure function lmoment_scalar_rank8(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function lmoment_scalar_rank8

    pure function lmoment_orders_rank8(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(8)
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        full_shape(1) = size(order)
        full_shape(2:8) = out_shape
        values = reshape(reduced, full_shape)
    end function lmoment_orders_rank8

    pure function lmoment_all_default_rank8(x, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, default_orders, sorted, standardize, &
            nan_policy, values)
    end function lmoment_all_default_rank8

    pure function lmoment_all_scalar_rank8(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        values = lmoment_all_scalar_assumed_rank(x, order, sorted, standardize, nan_policy)
    end function lmoment_all_scalar_rank8

    pure function lmoment_all_orders_rank8(x, order, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        call lmoment_all_orders_assumed_rank(x, order, sorted, standardize, nan_policy, values)
    end function lmoment_all_orders_rank8

    pure function lmoment_keepdims_default_rank8(x, sorted, standardize, nan_policy, axis) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: base_shape(8)
        integer :: full_shape(9)
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)
        integer :: d
        integer :: j

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_orders_assumed_rank(x, axis_value, default_orders, sorted, &
            standardize, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        base_shape(1) = size(default_orders)
        base_shape(2:8) = out_shape
        j = 1
        do d = 1, 9
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = base_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_default_rank8

    pure function lmoment_keepdims_scalar_rank8(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order !! positive L-moment order
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: full_shape(8)
        integer :: d
        integer :: j
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call lmoment_reduce_scalar_assumed_rank(x, axis_value, order, sorted, standardize, &
            nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        j = 1
        do d = 1, 8
            if (d == axis_value) then
                full_shape(d) = 1
            else
                full_shape(d) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_scalar_rank8

    pure function lmoment_keepdims_orders_rank8(x, order, axis, sorted, standardize, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! observations
        integer, intent(in) :: order(:) !! nonempty vector of positive L-moment orders
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        logical, intent(in), optional :: sorted !! true when each reduction slice is already sorted
        logical, intent(in), optional :: standardize !! standardize L-moments of order 3 and higher
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(9)
        integer :: j
        logical :: ok
        integer :: out_shape(7)
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
        j = 1
        do d = 1, 8
            if (d == axis_value) then
                full_shape(d + 1) = 1
            else
                full_shape(d + 1) = out_shape(j)
                j = j + 1
            end if
        end do
        values = reshape(reduced, full_shape)
    end function lmoment_keepdims_orders_rank8

end module scifort_descriptive_lmoment_1_8
