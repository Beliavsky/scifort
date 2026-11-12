! SPDX-License-Identifier: MIT
! Fixture driver for comparison with scipy.signal.argrelextrema.
program extrema_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal_relative_extrema, only : &
        argrelmax_mask_into, argrelmin_mask_into, argrelextrema_mask_into
    implicit none
    integer :: r, dims(15), selected, neighbors, mode_kind, op, code, n, ios, i
    real(dp), allocatable :: values(:)
    character(len=4) :: boundary

    do
        read(*,*,iostat=ios) r, dims, selected, neighbors, mode_kind, op
        if (ios < 0) exit
        if (ios /= 0) error stop 'unable to read extrema input header'
        if (r < 1 .or. r > 15 .or. any(dims(:r) < 0)) error stop 'invalid rank/dims'
        n = product(dims(:r))
        allocate(values(n))
        if (n > 0) read(*,*,iostat=ios) values
        if (ios > 0) error stop 'unable to read extrema signal'
        boundary = 'clip'
        if (mode_kind == 1) boundary = 'wrap'
        select case (r)
        case (1)
            block
                real(dp), allocatable :: x(:)
                logical, allocatable :: mask(:)
                allocate(x(dims(1)), mask(dims(1)))
                x = reshape(values, shape(x))
                select case (op)
                case (1)
                    call argrelmax_mask_into(x, mask, selected, neighbors, boundary, code)
                case (2)
                    call argrelmin_mask_into(x, mask, selected, neighbors, boundary, code)
                case (3)
                    call argrelextrema_mask_into(x, mask, greater_equal, &
                        selected, neighbors, boundary, code)
                end select
                write(*,'(i0,1x,i0)') code, count(mask)
                if (n > 0) write(*,'(*(i0,1x))') merge(1, 0, reshape(mask, [n]))
            end block
        case (2)
            block
                real(dp), allocatable :: x(:,:)
                logical, allocatable :: mask(:,:)
                allocate(x(dims(1),dims(2)), mask(dims(1),dims(2)))
                x = reshape(values, shape(x))
                select case (op)
                case (1)
                    call argrelmax_mask_into(x, mask, selected, neighbors, boundary, code)
                case (2)
                    call argrelmin_mask_into(x, mask, selected, neighbors, boundary, code)
                case (3)
                    call argrelextrema_mask_into(x, mask, greater_equal, &
                        selected, neighbors, boundary, code)
                end select
                write(*,'(i0,1x,i0)') code, count(mask)
                if (n > 0) write(*,'(*(i0,1x))') merge(1, 0, reshape(mask, [n]))
            end block
        case (3)
            block
                real(dp), allocatable :: x(:,:,:)
                logical, allocatable :: mask(:,:,:)
                allocate(x(dims(1),dims(2),dims(3)), mask(dims(1),dims(2),dims(3)))
                x = reshape(values, shape(x))
                select case (op)
                case (1)
                    call argrelmax_mask_into(x, mask, selected, neighbors, boundary, code)
                case (2)
                    call argrelmin_mask_into(x, mask, selected, neighbors, boundary, code)
                case (3)
                    call argrelextrema_mask_into(x, mask, greater_equal, &
                        selected, neighbors, boundary, code)
                end select
                write(*,'(i0,1x,i0)') code, count(mask)
                if (n > 0) write(*,'(*(i0,1x))') merge(1, 0, reshape(mask, [n]))
            end block
        case (4)
            block
                real(dp), allocatable :: x(:,:,:,:)
                logical, allocatable :: mask(:,:,:,:)
                allocate(x(dims(1),dims(2),dims(3),dims(4)), mask(dims(1),dims(2),dims(3),dims(4)))
                x = reshape(values, shape(x))
                select case (op)
                case (1)
                    call argrelmax_mask_into(x, mask, selected, neighbors, boundary, code)
                case (2)
                    call argrelmin_mask_into(x, mask, selected, neighbors, boundary, code)
                case (3)
                    call argrelextrema_mask_into(x, mask, greater_equal, &
                        selected, neighbors, boundary, code)
                end select
                write(*,'(i0,1x,i0)') code, count(mask)
                if (n > 0) write(*,'(*(i0,1x))') merge(1, 0, reshape(mask, [n]))
            end block
        case (5)
            block
                real(dp), allocatable :: x(:,:,:,:,:)
                logical, allocatable :: mask(:,:,:,:,:)
                allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5)), mask(dims(1),dims(2),dims(3),dims(4),dims(5)))
                x = reshape(values, shape(x))
                select case (op)
                case (1)
                    call argrelmax_mask_into(x, mask, selected, neighbors, boundary, code)
                case (2)
                    call argrelmin_mask_into(x, mask, selected, neighbors, boundary, code)
                case (3)
                    call argrelextrema_mask_into(x, mask, greater_equal, &
                        selected, neighbors, boundary, code)
                end select
                write(*,'(i0,1x,i0)') code, count(mask)
                if (n > 0) write(*,'(*(i0,1x))') merge(1, 0, reshape(mask, [n]))
            end block
        case (15)
            block
                real(dp), allocatable :: x(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
                logical, allocatable :: mask(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
                allocate(x(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8), &
                    dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15)))
                allocate(mask(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8), &
                    dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15)))
                x = reshape(values, shape(x))
                select case (op)
                case (1)
                    call argrelmax_mask_into(x, mask, selected, neighbors, boundary, code)
                case (2)
                    call argrelmin_mask_into(x, mask, selected, neighbors, boundary, code)
                case (3)
                    call argrelextrema_mask_into(x, mask, greater_equal, &
                        selected, neighbors, boundary, code)
                end select
                write(*,'(i0,1x,i0)') code, count(mask)
                if (n > 0) write(*,'(*(i0,1x))') merge(1, 0, reshape(mask, [n]))
            end block
        case default
            error stop 'unhandled test rank'
        end select
        deallocate(values)
    end do
contains
    pure logical function greater_equal(a, b)
        real(dp), intent(in) :: a !! candidate value
        real(dp), intent(in) :: b !! neighbor value
        greater_equal = a >= b
    end function greater_equal
end program extrema_nd_driver
