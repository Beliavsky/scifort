! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style relative extrema, with one-based Fortran indices.
! Independently implemented from the public argrelextrema specification.
module scifort_signal_relative_extrema
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: extrema_status_success = 0
    integer, parameter, public :: extrema_status_invalid_input = 1

    abstract interface
        pure logical function extrema_comparator(a, b)
            import dp
            real(dp), intent(in) :: a !! candidate value
            real(dp), intent(in) :: b !! neighbor value
        end function extrema_comparator
    end interface

    public :: argrelextrema, argrelmax, argrelmin
    public :: argrelextrema_mask_into, argrelmax_mask_into, argrelmin_mask_into
    public :: extrema_comparator

contains

    function argrelextrema(x, comparator, order, mode, status) result(indices)
        real(dp), intent(in) :: x(:) !! one-dimensional signal
        procedure(extrema_comparator) :: comparator !! pure scalar comparison predicate
        integer, intent(in), optional :: order !! number of neighbors on each side, >= 1
        character(*), intent(in), optional :: mode !! 'clip' (default) or 'wrap'
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        integer, allocatable :: indices(:) !! one-based extrema indices

        integer, allocatable :: scratch(:)
        integer :: i, j, n, nc, neighbor, distance
        logical :: keep, periodic

        allocate(indices(0))
        if (present(status)) status = extrema_status_invalid_input
        distance = 1
        if (present(order)) distance = order
        if (distance < 1) return
        periodic = .false.
        if (present(mode)) then
            select case (trim(mode))
            case ('clip')
            case ('wrap')
                periodic = .true.
            case default
                return
            end select
        end if

        n = size(x)
        nc = 0
        allocate(scratch(n))
        do i = 1, n
            keep = .true.
            do j = 1, distance
                if (periodic) then
                    ! Modulo avoids overflowing i+j or i-j for very large orders.
                    neighbor = mod( modulo(i - 1, n) + modulo(j, n), n) + 1
                else
                    neighbor = min(n, i + min(j, n))
                end if
                if (.not. comparator(x(i), x(neighbor))) then
                    keep = .false.
                    exit
                end if
                if (periodic) then
                    neighbor = modulo(modulo(i - 1, n) - modulo(j, n), n) + 1
                else
                    neighbor = max(1, i - min(j, n))
                end if
                if (.not. comparator(x(i), x(neighbor))) then
                    keep = .false.
                    exit
                end if
            end do
            if (.not. keep) cycle
            nc = nc + 1
            scratch(nc) = i
        end do
        indices = scratch(:nc)
        if (present(status)) status = extrema_status_success
    end function argrelextrema

    function argrelmax(x, order, mode, status) result(indices)
        real(dp), intent(in) :: x(:) !! one-dimensional signal
        integer, intent(in), optional :: order !! number of neighbors on each side
        character(*), intent(in), optional :: mode !! 'clip' or 'wrap'
        integer, intent(out), optional :: status !! zero on success
        integer, allocatable :: indices(:) !! one-based maxima indices

        indices = argrelextrema(x, strictly_greater, order, mode, status)
    end function argrelmax

    function argrelmin(x, order, mode, status) result(indices)
        real(dp), intent(in) :: x(:) !! one-dimensional signal
        integer, intent(in), optional :: order !! number of neighbors on each side
        character(*), intent(in), optional :: mode !! 'clip' or 'wrap'
        integer, intent(out), optional :: status !! zero on success
        integer, allocatable :: indices(:) !! one-based minima indices

        indices = argrelextrema(x, strictly_less, order, mode, status)
    end function argrelmin


    ! Logical masks are a SciFort extension to SciPy's coordinate-tuple output.
    ! This form supports all Fortran ranks, caller-owned storage, and strided arrays.
    subroutine argrelextrema_mask_into(x, mask, comparator, axis, order, mode, status)
        real(dp), intent(in) :: x(..) !! signal, rank one through fifteen
        logical, intent(inout) :: mask(..) !! output extrema mask, same shape as x
        procedure(extrema_comparator) :: comparator !! strict or user-defined comparison
        integer, intent(in), optional :: axis !! one-based axis, negative from the end, default 1
        integer, intent(in), optional :: order !! neighbors checked in each direction, default 1
        character(*), intent(in), optional :: mode !! endpoint behavior: 'clip' or 'wrap'
        integer, intent(out), optional :: status !! zero on success, nonzero on invalid arguments

        real(dp), allocatable :: values(:)
        logical, allocatable :: result_mask(:)
        integer, allocatable :: dims(:), strides(:)
        integer :: selected, distance, k, j, i, loc, idx, center, n, span
        character(len=:), allocatable :: boundary
        logical :: periodic, keep

        if (present(status)) status = extrema_status_invalid_input
        if (rank(x) < 1 .or. rank(x) > 15) return
        if (rank(mask) /= rank(x)) return
        if (any(shape(mask) /= shape(x))) return
        dims = shape(x)
        selected = 1
        if (present(axis)) selected = axis
        if (selected < 0) selected = selected + rank(x) + 1
        if (selected < 1 .or. selected > rank(x)) return
        distance = 1
        if (present(order)) distance = order
        if (distance < 1) return
        boundary = 'clip'
        if (present(mode)) boundary = trim(mode)
        select case (boundary)
        case ('clip')
            periodic = .false.
        case ('wrap')
            periodic = .true.
        case default
            return
        end select

        ! Reshape creates a contiguous copy of even a strided input section.
        ! Assignment back to mask happens only after all validation and work.
        select rank (x)
        rank (1)
            values = reshape(x, [size(x)])
        rank (2)
            values = reshape(x, [size(x)])
        rank (3)
            values = reshape(x, [size(x)])
        rank (4)
            values = reshape(x, [size(x)])
        rank (5)
            values = reshape(x, [size(x)])
        rank (6)
            values = reshape(x, [size(x)])
        rank (7)
            values = reshape(x, [size(x)])
        rank (8)
            values = reshape(x, [size(x)])
        rank (9)
            values = reshape(x, [size(x)])
        rank (10)
            values = reshape(x, [size(x)])
        rank (11)
            values = reshape(x, [size(x)])
        rank (12)
            values = reshape(x, [size(x)])
        rank (13)
            values = reshape(x, [size(x)])
        rank (14)
            values = reshape(x, [size(x)])
        rank (15)
            values = reshape(x, [size(x)])
        rank default
            return
        end select

        n = size(values)
        allocate(result_mask(n))
        result_mask = .false.
        if (n > 0) then
            span = dims(selected)
            if (span >= 1) then
                allocate(strides(size(dims)))
                strides(1) = 1
                do k = 2, size(dims)
                    strides(k) = strides(k - 1) * dims(k - 1)
                end do
                do i = 1, n
                    ! Coordinate of element i on the selected axis.
                    center = mod((i - 1) / strides(selected), span) + 1
                    keep = .true.
                    do j = 1, min(distance, span)
                        if (periodic) then
                            loc = mod(center - 1 + j, span) + 1
                        else
                            loc = min(span, center + j)
                        end if
                        idx = i + (loc - center) * strides(selected)
                        if (.not. comparator(values(i), values(idx))) then
                            keep = .false.
                            exit
                        end if
                        if (periodic) then
                            loc = modulo(center - 1 - j, span) + 1
                        else
                            loc = max(1, center - j)
                        end if
                        idx = i + (loc - center) * strides(selected)
                        if (.not. comparator(values(i), values(idx))) then
                            keep = .false.
                            exit
                        end if
                    end do
                    result_mask(i) = keep
                end do
            end if
        end if

        select rank (mask)
        rank (1)
            mask = reshape(result_mask, shape(mask))
        rank (2)
            mask = reshape(result_mask, shape(mask))
        rank (3)
            mask = reshape(result_mask, shape(mask))
        rank (4)
            mask = reshape(result_mask, shape(mask))
        rank (5)
            mask = reshape(result_mask, shape(mask))
        rank (6)
            mask = reshape(result_mask, shape(mask))
        rank (7)
            mask = reshape(result_mask, shape(mask))
        rank (8)
            mask = reshape(result_mask, shape(mask))
        rank (9)
            mask = reshape(result_mask, shape(mask))
        rank (10)
            mask = reshape(result_mask, shape(mask))
        rank (11)
            mask = reshape(result_mask, shape(mask))
        rank (12)
            mask = reshape(result_mask, shape(mask))
        rank (13)
            mask = reshape(result_mask, shape(mask))
        rank (14)
            mask = reshape(result_mask, shape(mask))
        rank (15)
            mask = reshape(result_mask, shape(mask))
        rank default
            return
        end select
        if (present(status)) status = extrema_status_success
    end subroutine argrelextrema_mask_into

    subroutine argrelmax_mask_into(x, mask, axis, order, mode, status)
        real(dp), intent(in) :: x(..) !! signal, rank one through fifteen
        logical, intent(inout) :: mask(..) !! location mask, with same shape as x
        integer, intent(in), optional :: axis !! one-based selected axis, default 1
        integer, intent(in), optional :: order !! strictly lower neighbors per side, default 1
        character(*), intent(in), optional :: mode !! 'clip' or 'wrap', default 'clip'
        integer, intent(out), optional :: status !! zero on success

        call argrelextrema_mask_into(x, mask, strictly_greater, axis, order, mode, status)
    end subroutine argrelmax_mask_into

    subroutine argrelmin_mask_into(x, mask, axis, order, mode, status)
        real(dp), intent(in) :: x(..) !! signal, rank one through fifteen
        logical, intent(inout) :: mask(..) !! location mask, with same shape as x
        integer, intent(in), optional :: axis !! one-based selected axis, default 1
        integer, intent(in), optional :: order !! strictly greater neighbors per side, default 1
        character(*), intent(in), optional :: mode !! 'clip' or 'wrap', default 'clip'
        integer, intent(out), optional :: status !! zero on success

        call argrelextrema_mask_into(x, mask, strictly_less, axis, order, mode, status)
    end subroutine argrelmin_mask_into

    pure logical function strictly_greater(a, b)
        real(dp), intent(in) :: a !! candidate value
        real(dp), intent(in) :: b !! neighbor value
        strictly_greater = a > b
    end function strictly_greater

    pure logical function strictly_less(a, b)
        real(dp), intent(in) :: a !! candidate value
        real(dp), intent(in) :: b !! neighbor value
        strictly_less = a < b
    end function strictly_less

end module scifort_signal_relative_extrema
