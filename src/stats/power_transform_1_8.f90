! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_power_transform_1_8
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    use scifort_power_transform_engine, only : boxcox_llf_all_assumed_rank, &
        boxcox_llf_reduce_assumed_rank, yeojohnson_llf_all_assumed_rank, &
        yeojohnson_llf_reduce_assumed_rank
    implicit none
    private

    interface boxcox_llf
        module procedure boxcox_llf_rank1
        module procedure boxcox_llf_rank2
        module procedure boxcox_llf_rank3
        module procedure boxcox_llf_rank4
        module procedure boxcox_llf_rank5
        module procedure boxcox_llf_rank6
        module procedure boxcox_llf_rank7
        module procedure boxcox_llf_rank8
    end interface boxcox_llf

    interface boxcox_llf_all
        module procedure boxcox_llf_all_rank1
        module procedure boxcox_llf_all_rank2
        module procedure boxcox_llf_all_rank3
        module procedure boxcox_llf_all_rank4
        module procedure boxcox_llf_all_rank5
        module procedure boxcox_llf_all_rank6
        module procedure boxcox_llf_all_rank7
        module procedure boxcox_llf_all_rank8
    end interface boxcox_llf_all

    interface boxcox_llf_keepdims
        module procedure boxcox_llf_keepdims_rank1
        module procedure boxcox_llf_keepdims_rank2
        module procedure boxcox_llf_keepdims_rank3
        module procedure boxcox_llf_keepdims_rank4
        module procedure boxcox_llf_keepdims_rank5
        module procedure boxcox_llf_keepdims_rank6
        module procedure boxcox_llf_keepdims_rank7
        module procedure boxcox_llf_keepdims_rank8
    end interface boxcox_llf_keepdims

    interface yeojohnson_llf
        module procedure yeojohnson_llf_rank1
        module procedure yeojohnson_llf_rank2
        module procedure yeojohnson_llf_rank3
        module procedure yeojohnson_llf_rank4
        module procedure yeojohnson_llf_rank5
        module procedure yeojohnson_llf_rank6
        module procedure yeojohnson_llf_rank7
        module procedure yeojohnson_llf_rank8
    end interface yeojohnson_llf

    interface yeojohnson_llf_all
        module procedure yeojohnson_llf_all_rank1
        module procedure yeojohnson_llf_all_rank2
        module procedure yeojohnson_llf_all_rank3
        module procedure yeojohnson_llf_all_rank4
        module procedure yeojohnson_llf_all_rank5
        module procedure yeojohnson_llf_all_rank6
        module procedure yeojohnson_llf_all_rank7
        module procedure yeojohnson_llf_all_rank8
    end interface yeojohnson_llf_all

    interface yeojohnson_llf_keepdims
        module procedure yeojohnson_llf_keepdims_rank1
        module procedure yeojohnson_llf_keepdims_rank2
        module procedure yeojohnson_llf_keepdims_rank3
        module procedure yeojohnson_llf_keepdims_rank4
        module procedure yeojohnson_llf_keepdims_rank5
        module procedure yeojohnson_llf_keepdims_rank6
        module procedure yeojohnson_llf_keepdims_rank7
        module procedure yeojohnson_llf_keepdims_rank8
    end interface yeojohnson_llf_keepdims

    public :: boxcox_llf, boxcox_llf_all, boxcox_llf_keepdims
    public :: yeojohnson_llf, yeojohnson_llf_all, yeojohnson_llf_keepdims

contains

    pure function boxcox_llf_rank1(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: values
        integer :: axis_value
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (ok .and. size(reduced) == 1) then
            values = reduced(1)
        else
            values = quiet_nan(0.0_dp)
        end if
    end function boxcox_llf_rank1

    pure function boxcox_llf_all_rank1(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank1

    pure function boxcox_llf_keepdims_rank1(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank1

    pure function boxcox_llf_rank2(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank2

    pure function boxcox_llf_all_rank2(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank2

    pure function boxcox_llf_keepdims_rank2(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank2

    pure function boxcox_llf_rank3(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank3

    pure function boxcox_llf_all_rank3(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank3

    pure function boxcox_llf_keepdims_rank3(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank3

    pure function boxcox_llf_rank4(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank4

    pure function boxcox_llf_all_rank4(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank4

    pure function boxcox_llf_keepdims_rank4(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank4

    pure function boxcox_llf_rank5(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank5

    pure function boxcox_llf_all_rank5(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank5

    pure function boxcox_llf_keepdims_rank5(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank5

    pure function boxcox_llf_rank6(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank6

    pure function boxcox_llf_all_rank6(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank6

    pure function boxcox_llf_keepdims_rank6(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank6

    pure function boxcox_llf_rank7(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank7

    pure function boxcox_llf_all_rank7(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank7

    pure function boxcox_llf_keepdims_rank7(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank7

    pure function boxcox_llf_rank8(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function boxcox_llf_rank8

    pure function boxcox_llf_all_rank8(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = boxcox_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function boxcox_llf_all_rank8

    pure function boxcox_llf_keepdims_rank8(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call boxcox_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function boxcox_llf_keepdims_rank8

    pure function yeojohnson_llf_rank1(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: values
        integer :: axis_value
        logical :: ok
        integer :: out_shape(0)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (ok .and. size(reduced) == 1) then
            values = reduced(1)
        else
            values = quiet_nan(0.0_dp)
        end if
    end function yeojohnson_llf_rank1

    pure function yeojohnson_llf_all_rank1(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank1

    pure function yeojohnson_llf_keepdims_rank1(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank1

    pure function yeojohnson_llf_rank2(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(1)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank2

    pure function yeojohnson_llf_all_rank2(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank2

    pure function yeojohnson_llf_keepdims_rank2(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank2

    pure function yeojohnson_llf_rank3(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(2)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank3

    pure function yeojohnson_llf_all_rank3(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank3

    pure function yeojohnson_llf_keepdims_rank3(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank3

    pure function yeojohnson_llf_rank4(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(3)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank4

    pure function yeojohnson_llf_all_rank4(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank4

    pure function yeojohnson_llf_keepdims_rank4(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank4

    pure function yeojohnson_llf_rank5(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(4)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank5

    pure function yeojohnson_llf_all_rank5(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank5

    pure function yeojohnson_llf_keepdims_rank5(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank5

    pure function yeojohnson_llf_rank6(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(5)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank6

    pure function yeojohnson_llf_all_rank6(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank6

    pure function yeojohnson_llf_keepdims_rank6(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank6

    pure function yeojohnson_llf_rank7(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(6)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank7

    pure function yeojohnson_llf_all_rank7(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank7

    pure function yeojohnson_llf_keepdims_rank7(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank7

    pure function yeojohnson_llf_rank8(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
        real(dp), allocatable :: values(:,:,:,:,:,:,:)
        integer :: axis_value
        logical :: ok
        integer :: out_shape(7)
        real(dp), allocatable :: reduced(:)
        axis_value = 1
        if (present(axis)) axis_value = axis
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
        if (.not. ok) then
            allocate(values(0,0,0,0,0,0,0))
            return
        end if
        values = reshape(reduced, out_shape)
    end function yeojohnson_llf_rank8

    pure function yeojohnson_llf_all_rank8(lmbda, data, nan_policy) result(value)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:)
        character(len=*), intent(in), optional :: nan_policy
        real(dp) :: value
        value = yeojohnson_llf_all_assumed_rank(lmbda, data, nan_policy)
    end function yeojohnson_llf_all_rank8

    pure function yeojohnson_llf_keepdims_rank8(lmbda, data, axis, nan_policy) result(values)
        real(dp), intent(in) :: lmbda
        real(dp), intent(in) :: data(:,:,:,:,:,:,:,:)
        integer, intent(in), optional :: axis
        character(len=*), intent(in), optional :: nan_policy
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
        call yeojohnson_llf_reduce_assumed_rank(lmbda, data, axis_value, nan_policy, reduced, out_shape, ok)
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
    end function yeojohnson_llf_keepdims_rank8

end module scifort_power_transform_1_8
