! SPDX-License-Identifier: MIT
! Numerical reference driver: input family, length, symmetric, p1, p2, norm.
program window_advanced_driver
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window, chebwin, taylor, kaiser_bessel_derived
    implicit none
    real(dp), allocatable :: w(:)
    real(dp) :: p1, p2
    integer :: family, n, sy, normal, st, ios
    do
        read(*,*,iostat=ios) family, n, sy, p1, p2, normal
        if (ios /= 0) exit
        select case (family)
        case (1)
            w = chebwin(n, p1, sym=sy==1, status=st)
        case (2)
            w = taylor(n, nbar=int(p1), sll=p2, norm=normal==1, sym=sy==1, status=st)
        case (3)
            w = kaiser_bessel_derived(n, p1, sym=sy==1, status=st)
        case (4)
            w = get_window('chebwin', n, parameter=p1, fftbins=sy/=1, status=st)
        case (5)
            w = get_window('taylor', n, parameter=p1, parameter2=p2, fftbins=sy/=1, status=st)
        case (6)
            w = get_window('kbd', n, parameter=p1, fftbins=sy/=1, status=st)
        case default
            error stop 'invalid family'
        end select
        write(*,'(i0,1x,i0,*(1x,es25.17))') st, size(w), w
    end do
end program window_advanced_driver
