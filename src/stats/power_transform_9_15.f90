! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_power_transform_9_15
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_power_transform_engine, only : boxcox_llf_all_assumed_rank, &
        boxcox_llf_reduce_assumed_rank, yeojohnson_llf_all_assumed_rank, &
        yeojohnson_llf_reduce_assumed_rank
    implicit none
    private

    interface boxcox_llf
        module procedure boxcox_llf_rank9
        module procedure boxcox_llf_rank10
        module procedure boxcox_llf_rank11
        module procedure boxcox_llf_rank12
        module procedure boxcox_llf_rank13
        module procedure boxcox_llf_rank14
        module procedure boxcox_llf_rank15
    end interface boxcox_llf

    interface boxcox_llf_all
        module procedure boxcox_llf_all_rank9
        module procedure boxcox_llf_all_rank10
        module procedure boxcox_llf_all_rank11
        module procedure boxcox_llf_all_rank12
        module procedure boxcox_llf_all_rank13
        module procedure boxcox_llf_all_rank14
        module procedure boxcox_llf_all_rank15
    end interface boxcox_llf_all

    interface boxcox_llf_keepdims
        module procedure boxcox_llf_keepdims_rank9
        module procedure boxcox_llf_keepdims_rank10
        module procedure boxcox_llf_keepdims_rank11
        module procedure boxcox_llf_keepdims_rank12
        module procedure boxcox_llf_keepdims_rank13
        module procedure boxcox_llf_keepdims_rank14
        module procedure boxcox_llf_keepdims_rank15
    end interface boxcox_llf_keepdims

    interface yeojohnson_llf
        module procedure yeojohnson_llf_rank9
        module procedure yeojohnson_llf_rank10
        module procedure yeojohnson_llf_rank11
        module procedure yeojohnson_llf_rank12
        module procedure yeojohnson_llf_rank13
        module procedure yeojohnson_llf_rank14
        module procedure yeojohnson_llf_rank15
    end interface yeojohnson_llf

    interface yeojohnson_llf_all
        module procedure yeojohnson_llf_all_rank9
        module procedure yeojohnson_llf_all_rank10
        module procedure yeojohnson_llf_all_rank11
        module procedure yeojohnson_llf_all_rank12
        module procedure yeojohnson_llf_all_rank13
        module procedure yeojohnson_llf_all_rank14
        module procedure yeojohnson_llf_all_rank15
    end interface yeojohnson_llf_all

    interface yeojohnson_llf_keepdims
        module procedure yeojohnson_llf_keepdims_rank9
        module procedure yeojohnson_llf_keepdims_rank10
        module procedure yeojohnson_llf_keepdims_rank11
        module procedure yeojohnson_llf_keepdims_rank12
        module procedure yeojohnson_llf_keepdims_rank13
        module procedure yeojohnson_llf_keepdims_rank14
        module procedure yeojohnson_llf_keepdims_rank15
    end interface yeojohnson_llf_keepdims

    public :: boxcox_llf, boxcox_llf_all, boxcox_llf_keepdims
    public :: yeojohnson_llf, yeojohnson_llf_all, yeojohnson_llf_keepdims

contains

    pure function boxcox_llf_rank9(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank9

    pure function boxcox_llf_all_rank9(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank9

    pure function boxcox_llf_keepdims_rank9(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank9

    pure function boxcox_llf_rank10(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank10

    pure function boxcox_llf_all_rank10(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank10

    pure function boxcox_llf_keepdims_rank10(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank10

    pure function boxcox_llf_rank11(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank11

    pure function boxcox_llf_all_rank11(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank11

    pure function boxcox_llf_keepdims_rank11(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank11

    pure function boxcox_llf_rank12(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank12

    pure function boxcox_llf_all_rank12(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank12

    pure function boxcox_llf_keepdims_rank12(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank12

    pure function boxcox_llf_rank13(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank13

    pure function boxcox_llf_all_rank13(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank13

    pure function boxcox_llf_keepdims_rank13(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank13

    pure function boxcox_llf_rank14(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank14

    pure function boxcox_llf_all_rank14(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank14

    pure function boxcox_llf_keepdims_rank14(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank14

    pure function boxcox_llf_rank15(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank15

    pure function boxcox_llf_all_rank15(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank15

    pure function boxcox_llf_keepdims_rank15(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank15

    pure function yeojohnson_llf_rank9(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(8)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank9

    pure function yeojohnson_llf_all_rank9(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank9

    pure function yeojohnson_llf_keepdims_rank9(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank9

    pure function yeojohnson_llf_rank10(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(9)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank10

    pure function yeojohnson_llf_all_rank10(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank10

    pure function yeojohnson_llf_keepdims_rank10(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank10

    pure function yeojohnson_llf_rank11(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(10)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank11

    pure function yeojohnson_llf_all_rank11(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank11

    pure function yeojohnson_llf_keepdims_rank11(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank11

    pure function yeojohnson_llf_rank12(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(11)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank12

    pure function yeojohnson_llf_all_rank12(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank12

    pure function yeojohnson_llf_keepdims_rank12(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank12

    pure function yeojohnson_llf_rank13(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(12)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank13

    pure function yeojohnson_llf_all_rank13(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank13

    pure function yeojohnson_llf_keepdims_rank13(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank13

    pure function yeojohnson_llf_rank14(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(13)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank14

    pure function yeojohnson_llf_all_rank14(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank14

    pure function yeojohnson_llf_keepdims_rank14(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank14

    pure function yeojohnson_llf_rank15(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(14)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0,0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank15

    pure function yeojohnson_llf_all_rank15(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank15

    pure function yeojohnson_llf_keepdims_rank15(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank15

end module scifort_power_transform_9_15
