! SPDX-License-Identifier: MIT
! Text protocol used by the independent SciPy 1.17 numerical comparison.
program kaiser_spline_driver
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : gauss_spline
    use scifort_signal_design, only : kaiser_beta, kaiser_atten, kaiserord, &
        kaiser_order_result, firwin
    implicit none
    character(1) :: op
    real(dp) :: x, y
    real(dp), allocatable :: taps(:)
    integer :: ios, n, st
    type(kaiser_order_result) :: spec

    do
        read(*,*,iostat=ios) op, x, y, n
        if (ios /= 0) exit
        select case(op)
        case('g')
            print '(es26.17e3)', gauss_spline(x,n)
        case('b')
            print '(es26.17e3)', kaiser_beta(x)
        case('a')
            print '(es26.17e3)', kaiser_atten(n,x)
        case('o')
            spec = kaiserord(x,y)
            print '(i0,1x,i0,1x,es26.17e3)', spec%status, spec%numtaps, spec%beta
        case('f')
            taps = firwin(n,[x],window='kaiser',window_parameter=y,status=st)
            print '(i0,1x,i0)',st,size(taps)
            print '(*(1x,es26.17e3))',taps
        end select
    end do
end program kaiser_spline_driver
