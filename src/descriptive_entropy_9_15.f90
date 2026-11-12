! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_descriptive_entropy_9_15
    use scifort_descriptive_entropy_engine, only : differential_entropy_all_assumed_rank, &
        differential_entropy_reduce_assumed_rank, entropy_all_assumed_rank, &
        entropy_reduce_assumed_rank
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    interface entropy
        module procedure entropy_rank9
        module procedure entropy_rank10
        module procedure entropy_rank11
        module procedure entropy_rank12
        module procedure entropy_rank13
        module procedure entropy_rank14
        module procedure entropy_rank15
    end interface entropy

    interface entropy_all
        module procedure entropy_all_rank9
        module procedure entropy_all_rank10
        module procedure entropy_all_rank11
        module procedure entropy_all_rank12
        module procedure entropy_all_rank13
        module procedure entropy_all_rank14
        module procedure entropy_all_rank15
    end interface entropy_all

    interface entropy_keepdims
        module procedure entropy_keepdims_rank9
        module procedure entropy_keepdims_rank10
        module procedure entropy_keepdims_rank11
        module procedure entropy_keepdims_rank12
        module procedure entropy_keepdims_rank13
        module procedure entropy_keepdims_rank14
        module procedure entropy_keepdims_rank15
    end interface entropy_keepdims

    interface differential_entropy
        module procedure differential_entropy_rank9
        module procedure differential_entropy_rank10
        module procedure differential_entropy_rank11
        module procedure differential_entropy_rank12
        module procedure differential_entropy_rank13
        module procedure differential_entropy_rank14
        module procedure differential_entropy_rank15
    end interface differential_entropy

    interface differential_entropy_all
        module procedure differential_entropy_all_rank9
        module procedure differential_entropy_all_rank10
        module procedure differential_entropy_all_rank11
        module procedure differential_entropy_all_rank12
        module procedure differential_entropy_all_rank13
        module procedure differential_entropy_all_rank14
        module procedure differential_entropy_all_rank15
    end interface differential_entropy_all

    interface differential_entropy_keepdims
        module procedure differential_entropy_keepdims_rank9
        module procedure differential_entropy_keepdims_rank10
        module procedure differential_entropy_keepdims_rank11
        module procedure differential_entropy_keepdims_rank12
        module procedure differential_entropy_keepdims_rank13
        module procedure differential_entropy_keepdims_rank14
        module procedure differential_entropy_keepdims_rank15
    end interface differential_entropy_keepdims

    public :: entropy
    public :: entropy_all
    public :: entropy_keepdims
    public :: differential_entropy
    public :: differential_entropy_all
    public :: differential_entropy_keepdims

contains

    pure function entropy_rank9(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0))
        end if
    end function entropy_rank9

    pure function entropy_all_rank9(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank9

    pure function entropy_keepdims_rank9(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(9)
        logical :: ok
        integer :: out_shape(8)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 9
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank9

    pure function differential_entropy_rank9(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0))
        end if
    end function differential_entropy_rank9

    pure function differential_entropy_all_rank9(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank9

    pure function differential_entropy_keepdims_rank9(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(9)
        logical :: ok
        integer :: out_shape(8)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 9
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank9

    pure function entropy_rank10(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0))
        end if
    end function entropy_rank10

    pure function entropy_all_rank10(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank10

    pure function entropy_keepdims_rank10(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(10)
        logical :: ok
        integer :: out_shape(9)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 10
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank10

    pure function differential_entropy_rank10(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0))
        end if
    end function differential_entropy_rank10

    pure function differential_entropy_all_rank10(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank10

    pure function differential_entropy_keepdims_rank10(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(10)
        logical :: ok
        integer :: out_shape(9)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 10
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank10

    pure function entropy_rank11(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0))
        end if
    end function entropy_rank11

    pure function entropy_all_rank11(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank11

    pure function entropy_keepdims_rank11(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(11)
        logical :: ok
        integer :: out_shape(10)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 11
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank11

    pure function differential_entropy_rank11(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0))
        end if
    end function differential_entropy_rank11

    pure function differential_entropy_all_rank11(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank11

    pure function differential_entropy_keepdims_rank11(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(11)
        logical :: ok
        integer :: out_shape(10)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 11
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank11

    pure function entropy_rank12(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
        end if
    end function entropy_rank12

    pure function entropy_all_rank12(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank12

    pure function entropy_keepdims_rank12(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(12)
        logical :: ok
        integer :: out_shape(11)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 12
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank12

    pure function differential_entropy_rank12(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
        end if
    end function differential_entropy_rank12

    pure function differential_entropy_all_rank12(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank12

    pure function differential_entropy_keepdims_rank12(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(12)
        logical :: ok
        integer :: out_shape(11)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 12
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank12

    pure function entropy_rank13(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
        end if
    end function entropy_rank13

    pure function entropy_all_rank13(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank13

    pure function entropy_keepdims_rank13(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(13)
        logical :: ok
        integer :: out_shape(12)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 13
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank13

    pure function differential_entropy_rank13(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
        end if
    end function differential_entropy_rank13

    pure function differential_entropy_all_rank13(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank13

    pure function differential_entropy_keepdims_rank13(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(13)
        logical :: ok
        integer :: out_shape(12)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 13
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank13

    pure function entropy_rank14(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
        end if
    end function entropy_rank14

    pure function entropy_all_rank14(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank14

    pure function entropy_keepdims_rank14(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(14)
        logical :: ok
        integer :: out_shape(13)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 14
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank14

    pure function differential_entropy_rank14(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
        end if
    end function differential_entropy_rank14

    pure function differential_entropy_all_rank14(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank14

    pure function differential_entropy_keepdims_rank14(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(14)
        logical :: ok
        integer :: out_shape(13)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 14
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank14

    pure function entropy_rank15(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        end if
    end function entropy_rank15

    pure function entropy_all_rank15(pk, qk, base, nan_policy) result(value)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = entropy_all_assumed_rank(pk, qk, base, nan_policy)
    end function entropy_all_rank15

    pure function entropy_keepdims_rank15(pk, qk, base, axis, nan_policy) result(values)
        real(dp), intent(in) :: pk(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! probability weights
        real(dp), intent(in), optional :: qk(..) !! broadcastable comparison probabilities
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(15)
        logical :: ok
        integer :: out_shape(14)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call entropy_reduce_assumed_rank(pk, qk, axis_value, base, nan_policy, &
            reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 15
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function entropy_keepdims_rank15

    pure function differential_entropy_rank15(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        logical :: ok
        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (ok) then
            values = reshape(reduced, out_shape)
        else
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
        end if
    end function differential_entropy_rank15

    pure function differential_entropy_all_rank15(x, window_length, base, method, nan_policy) result(value)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp) :: value

        value = differential_entropy_all_assumed_rank(x, window_length, base, method, nan_policy)
    end function differential_entropy_all_rank15

    pure function differential_entropy_keepdims_rank15(x, window_length, base, axis, method, nan_policy) result(values)
        real(dp), intent(in) :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! continuous sample values
        integer, intent(in), optional :: window_length !! spacing window length
        real(dp), intent(in), optional :: base !! logarithm base; natural when absent
        integer, intent(in), optional :: axis !! one-based Fortran reduction axis, default 1
        character(len=*), intent(in), optional :: method !! vasicek/van es/ebrahimi/correa/auto
        character(len=*), intent(in), optional :: nan_policy !! propagate/omit/raise
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)

        integer :: axis_value
        integer :: d
        integer :: full_shape(15)
        logical :: ok
        integer :: out_shape(14)
        integer :: out_dim
        real(dp), allocatable :: reduced(:)

        axis_value = 1
        if (present(axis)) axis_value = axis
        call differential_entropy_reduce_assumed_rank(x, axis_value, window_length, &
            base, method, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        out_dim = 0
        do d = 1, 15
            if (d == axis_value) then
                full_shape(d) = 1
            else
                out_dim = out_dim + 1
                full_shape(d) = out_shape(out_dim)
            end if
        end do
        values = reshape(reduced, full_shape)
    end function differential_entropy_keepdims_rank15

end module scifort_descriptive_entropy_9_15
