! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Stream driver for independently comparing SciFort and SciPy extrema results.
program extrema_excursion_driver
    use scifort_kinds, only : dp
    use scifort_signal_relative_extrema, only : argrelmax, argrelmin, argrelextrema
    use scifort_signal_peak_excursion, only : peak_detect
    implicit none

    integer :: n, order, periodic, io
    character(4) :: mode
    real(dp) :: delta
    real(dp), allocatable :: x(:), mx(:), mn(:)
    integer, allocatable :: imx(:), imn(:), ind(:)

    do
        read(*, *, iostat=io) n, order, periodic, delta
        if (io /= 0) exit
        if (n < 0) error stop 'invalid signal length'
        allocate(x(n))
        if (n > 0) read(*, *) x
        mode = 'clip'
        if (periodic /= 0) mode = 'wrap'
        ind = argrelmax(x, order, mode)
        call write_indices(ind)
        ind = argrelmin(x, order, mode)
        call write_indices(ind)
        ind = argrelextrema(x, greater_or_equal, order, mode)
        call write_indices(ind)
        call peak_detect(x, delta, imx, mx, imn, mn)
        call write_indices(imx)
        call write_values(mx)
        call write_indices(imn)
        call write_values(mn)
        deallocate(x)
    end do

contains

    pure logical function greater_or_equal(a, b)
        real(dp), intent(in) :: a !! candidate
        real(dp), intent(in) :: b !! neighbor
        greater_or_equal = a >= b
    end function greater_or_equal

    subroutine write_indices(a)
        integer, intent(in) :: a(:) !! indices
        print '(i0,*(1x,i0))', size(a), a
    end subroutine write_indices

    subroutine write_values(a)
        real(dp), intent(in) :: a(:) !! values
        print '(i0,*(1x,es26.17e3))', size(a), a
    end subroutine write_values

end program extrema_excursion_driver
