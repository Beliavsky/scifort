! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program czt_driver
    use scifort_kinds, only : dp
    use scifort_signal_czt, only : czt, zoom_fft
    implicit none
    integer :: count, k, n, m, mode, endpt, st
    real(dp) :: fn1, fn2, fs, ar, ai, wr, wi
    real(dp), allocatable :: xr(:), xi(:)
    complex(dp), allocatable :: y(:)

    read(*,*) count
    do k = 1,count
        read(*,*) n, m, mode
        read(*,*) fn1, fn2, fs, endpt
        read(*,*) ar, ai, wr, wi
        allocate(xr(n),xi(n))
        read(*,*) xr
        read(*,*) xi
        select case(mode)
        case(1)
            y = czt(cmplx(xr,xi,dp),m=m,status=st)
        case(2)
            y = czt(cmplx(xr,xi,dp),m=m,a=cmplx(ar,ai,dp),w=cmplx(wr,wi,dp),status=st)
        case(3)
            y = zoom_fft(cmplx(xr,xi,dp),[fn1,fn2],m=m,fs=fs,endpoint=endpt==1,status=st)
        case(4)
            y = zoom_fft(cmplx(xr,xi,dp),fn2,m=m,fs=fs,endpoint=endpt==1,status=st)
        end select
        print '(2(i0,1x))', st, size(y)
        if (size(y) > 0) then
            print '(*(es25.16,1x))', real(y,dp)
            print '(*(es25.16,1x))', aimag(y)
        end if
        deallocate(xr,xi,y)
    end do
end program czt_driver
