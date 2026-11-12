! SPDX-License-Identifier: MIT
! SciPy reference driver for the independent polyphase filtering implementation.
program resample_driver
    use scifort_signal_resample, only : upfirdn, resample_poly, decimate
    use scifort_kinds, only : dp
    implicit none
    integer :: n, m, nu, nd, stat, kind, i, j, order, centered, ierr
    character(32) :: mode, window
    real(dp), allocatable :: a(:), b(:), v(:), h(:)
    complex(dp), allocatable :: z(:), cz(:), hc(:)
    do
        read (*,*,iostat=ierr) kind,n,m,nu,nd,mode,window,order,centered
        if (ierr/=0) exit
        allocate(a(n),b(n),v(0),h(m),z(n),hc(m),cz(0))
        read (*,*) (a(i),i=1,n)
        read (*,*) (b(i),i=1,n)
        if (m>0) read (*,*) (h(i),i=1,m)
        z=cmplx(a,b,dp)
        if (m>0) hc=cmplx(h,0.0_dp,dp)
        select case(kind)
        case(1)
            v=upfirdn(h,a,up=nu,down=nd,mode=mode,status=stat)
        case(2)
            cz=upfirdn(h,z,up=nu,down=nd,mode=mode,status=stat)
        case(3)
            cz=upfirdn(hc,z,up=nu,down=nd,mode=mode,status=stat)
        case(4)
            if (m>0) then
                v=resample_poly(a,nu,nd,taps=h,padtype=mode,status=stat)
            else
                v=resample_poly(a,nu,nd,window=window,padtype=mode,status=stat)
            end if
        case(5)
            if (m>0) then
                cz=resample_poly(z,nu,nd,taps=h,padtype=mode,status=stat)
            else
                cz=resample_poly(z,nu,nd,window=window,padtype=mode,status=stat)
            end if
        case(6)
            v=decimate(a,nd,n=order,ftype='fir',zero_phase=centered==1,status=stat)
        case(7)
            cz=decimate(z,nd,n=order,ftype='fir',zero_phase=centered==1,status=stat)
        end select
        if (kind==1 .or. kind==4 .or. kind==6) then
            print '(i0,1x,i0)', stat,size(v)
            do i=1,size(v)
                print '(es25.17e3)', v(i)
            end do
        else
            print '(i0,1x,i0)', stat,size(cz)
            do i=1,size(cz)
                print '(2(es25.17e3,1x))',real(cz(i),dp),aimag(cz(i))
            end do
        end if
        deallocate(a,b,v,h,z,hc,cz)
    end do
end program resample_driver
