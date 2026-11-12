! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_entropy_1_8
    use scifort_descriptive_entropy_engine, only : differential_entropy_all_assumed_rank, &
        differential_entropy_reduce_assumed_rank, entropy_all_assumed_rank, &
        entropy_reduce_assumed_rank
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface entropy
        module procedure entropy_rank1
        module procedure entropy_rank2
        module procedure entropy_rank3
        module procedure entropy_rank4
        module procedure entropy_rank5
        module procedure entropy_rank6
        module procedure entropy_rank7
        module procedure entropy_rank8
    end interface entropy

    interface entropy_all
        module procedure entropy_all_rank1
        module procedure entropy_all_rank2
        module procedure entropy_all_rank3
        module procedure entropy_all_rank4
        module procedure entropy_all_rank5
        module procedure entropy_all_rank6
        module procedure entropy_all_rank7
        module procedure entropy_all_rank8
    end interface entropy_all

    interface entropy_keepdims
        module procedure entropy_keepdims_rank1
        module procedure entropy_keepdims_rank2
        module procedure entropy_keepdims_rank3
        module procedure entropy_keepdims_rank4
        module procedure entropy_keepdims_rank5
        module procedure entropy_keepdims_rank6
        module procedure entropy_keepdims_rank7
        module procedure entropy_keepdims_rank8
    end interface entropy_keepdims

    interface differential_entropy
        module procedure differential_entropy_rank1
        module procedure differential_entropy_rank2
        module procedure differential_entropy_rank3
        module procedure differential_entropy_rank4
        module procedure differential_entropy_rank5
        module procedure differential_entropy_rank6
        module procedure differential_entropy_rank7
        module procedure differential_entropy_rank8
    end interface differential_entropy

    interface differential_entropy_all
        module procedure differential_entropy_all_rank1
        module procedure differential_entropy_all_rank2
        module procedure differential_entropy_all_rank3
        module procedure differential_entropy_all_rank4
        module procedure differential_entropy_all_rank5
        module procedure differential_entropy_all_rank6
        module procedure differential_entropy_all_rank7
        module procedure differential_entropy_all_rank8
    end interface differential_entropy_all

    interface differential_entropy_keepdims
        module procedure differential_entropy_keepdims_rank1
        module procedure differential_entropy_keepdims_rank2
        module procedure differential_entropy_keepdims_rank3
        module procedure differential_entropy_keepdims_rank4
        module procedure differential_entropy_keepdims_rank5
        module procedure differential_entropy_keepdims_rank6
        module procedure differential_entropy_keepdims_rank7
        module procedure differential_entropy_keepdims_rank8
    end interface differential_entropy_keepdims

    public :: entropy
    public :: entropy_all
    public :: entropy_keepdims
    public :: differential_entropy
    public :: differential_entropy_all
    public :: differential_entropy_keepdims

contains

    pure function entropy_rank1(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        integer :: axis_value
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok .and. size(reduced) == 1) then
            values = reduced(1)
        else
            values = quiet_nan(0.0_dp)
        end if
    end function entropy_rank1

    pure function entropy_all_rank1(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank1

    pure function entropy_keepdims_rank1(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(1)
        logical :: ok
        integer :: out_shape(0)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        out_dim = 0
        do d = 1, 1
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank1

    pure function differential_entropy_rank1(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: values

        integer :: axis_value
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok .and. size(reduced) == 1) then
            values = reduced(1)
        else
            values = quiet_nan(0.0_dp)
        end if
    end function differential_entropy_rank1

    pure function differential_entropy_all_rank1(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank1

    pure function differential_entropy_keepdims_rank1(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(1)
        logical :: ok
        integer :: out_shape(0)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        out_dim = 0
        do d = 1, 1
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank1

    pure function entropy_rank2(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0))
        end if
    end function entropy_rank2

    pure function entropy_all_rank2(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank2

    pure function entropy_keepdims_rank2(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(2)
        logical :: ok
        integer :: out_shape(1)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        out_dim = 0
        do d = 1, 2
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank2

    pure function differential_entropy_rank2(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0))
        end if
    end function differential_entropy_rank2

    pure function differential_entropy_all_rank2(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank2

    pure function differential_entropy_keepdims_rank2(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(2)
        logical :: ok
        integer :: out_shape(1)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        out_dim = 0
        do d = 1, 2
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank2

    pure function entropy_rank3(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0))
        end if
    end function entropy_rank3

    pure function entropy_all_rank3(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank3

    pure function entropy_keepdims_rank3(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(3)
        logical :: ok
        integer :: out_shape(2)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 3
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank3

    pure function differential_entropy_rank3(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0))
        end if
    end function differential_entropy_rank3

    pure function differential_entropy_all_rank3(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank3

    pure function differential_entropy_keepdims_rank3(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(3)
        logical :: ok
        integer :: out_shape(2)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 3
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank3

    pure function entropy_rank4(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0))
        end if
    end function entropy_rank4

    pure function entropy_all_rank4(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank4

    pure function entropy_keepdims_rank4(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(4)
        logical :: ok
        integer :: out_shape(3)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 4
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank4

    pure function differential_entropy_rank4(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0))
        end if
    end function differential_entropy_rank4

    pure function differential_entropy_all_rank4(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank4

    pure function differential_entropy_keepdims_rank4(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(4)
        logical :: ok
        integer :: out_shape(3)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 4
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank4

    pure function entropy_rank5(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0))
        end if
    end function entropy_rank5

    pure function entropy_all_rank5(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank5

    pure function entropy_keepdims_rank5(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(5)
        logical :: ok
        integer :: out_shape(4)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 5
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank5

    pure function differential_entropy_rank5(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0))
        end if
    end function differential_entropy_rank5

    pure function differential_entropy_all_rank5(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank5

    pure function differential_entropy_keepdims_rank5(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(5)
        logical :: ok
        integer :: out_shape(4)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 5
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank5

    pure function entropy_rank6(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0))
        end if
    end function entropy_rank6

    pure function entropy_all_rank6(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank6

    pure function entropy_keepdims_rank6(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(6)
        logical :: ok
        integer :: out_shape(5)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 6
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank6

    pure function differential_entropy_rank6(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0))
        end if
    end function differential_entropy_rank6

    pure function differential_entropy_all_rank6(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank6

    pure function differential_entropy_keepdims_rank6(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(6)
        logical :: ok
        integer :: out_shape(5)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 6
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank6

    pure function entropy_rank7(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0))
        end if
    end function entropy_rank7

    pure function entropy_all_rank7(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank7

    pure function entropy_keepdims_rank7(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(7)
        logical :: ok
        integer :: out_shape(6)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 7
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank7

    pure function differential_entropy_rank7(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0))
        end if
    end function differential_entropy_rank7

    pure function differential_entropy_all_rank7(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank7

    pure function differential_entropy_keepdims_rank7(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(7)
        logical :: ok
        integer :: out_shape(6)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 7
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank7

    pure function entropy_rank8(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0))
        end if
    end function entropy_rank8

    pure function entropy_all_rank8(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank8

    pure function entropy_keepdims_rank8(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(8)
        logical :: ok
        integer :: out_shape(7)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 8
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank8

    pure function differential_entropy_rank8(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0))
        end if
    end function differential_entropy_rank8

    pure function differential_entropy_all_rank8(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank8

    pure function differential_entropy_keepdims_rank8(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(8)
        logical :: ok
        integer :: out_shape(7)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 8
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank8

end module scifort_descriptive_entropy_1_8
