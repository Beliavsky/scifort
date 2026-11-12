! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Standalone comparison driver. Compile with scifort_kinds,
! scifort_signal_windows_extra, and scifort_signal_spectral.
program window_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window
    implicit none
    character(64) :: name
    real(dp) :: p1, p2
    real(dp), allocatable :: w(:)
    integer :: j, nc, n, ip, h1, h2, st, i
    logical :: periodic

    read(*,*) nc
    do j=1,nc
        read(*,*) name, n, ip, h1, p1, h2, p2
        periodic = ip /= 0
        if (h1 /= 0) then
            if (h2 /= 0) then
                w = get_window(trim(name), n, fftbins=periodic, &
                               status=st, parameter=p1, parameter2=p2)
            else
                w = get_window(trim(name), n, fftbins=periodic, &
                               status=st, parameter=p1)
            end if
        else
            w = get_window(trim(name), n, fftbins=periodic, status=st)
        end if
        print '(2(i0,1x))', st, size(w)
        do i=1,size(w)
            print '(es26.17e3)', w(i)
        end do
    end do
end program window_reference_driver
