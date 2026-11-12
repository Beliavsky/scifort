! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Test-only streaming driver for verify_filter_diagnostics.py.
program filter_diagnostics_driver
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfiltic
    use scifort_signal_design, only : group_delay, group_delay_result
    implicit none

    real(dp), allocatable :: b(:), a(:), y(:), x(:), w(:), zi(:)
    complex(dp), allocatable :: bc(:), ac(:), yc(:), xc(:), zic(:)
    type(group_delay_result) :: gd
    integer :: operation, nb, na, ny, nx, nw, status, ios, i
    real(dp) :: fs

    do
        read(*, *, iostat=ios) operation, nb, na, ny, nx, nw, fs
        if (ios /= 0) exit
        select case (operation)
        case (1,3)
            allocate(b(nb),a(na))
            read(*,*) b
            read(*,*) a
        case (2,4)
            allocate(bc(nb),ac(na))
            read(*,*) bc
            read(*,*) ac
        end select
        select case (operation)
        case (1)
            allocate(y(ny))
            if (ny > 0) read(*,*) y
            if (nx >= 0) then
                allocate(x(nx))
                if (nx > 0) read(*,*) x
                zi = lfiltic(b,a,y,x,status)
                deallocate(x)
            else
                zi = lfiltic(b,a,y,status=status)
            end if
            print '(i0,1x,i0)',status,size(zi)
            if (size(zi)>0) print '(*(es26.17e3,1x))',zi
            deallocate(b,a,y,zi)
        case (2)
            allocate(yc(ny))
            if (ny > 0) read(*,*) yc
            if (nx >= 0) then
                allocate(xc(nx))
                if (nx > 0) read(*,*) xc
                zic = lfiltic(bc,ac,yc,xc,status)
                deallocate(xc)
            else
                zic = lfiltic(bc,ac,yc,status=status)
            end if
            print '(i0,1x,i0)',status,size(zic)
            do i=1,size(zic)
                print '(2(es26.17e3,1x))',real(zic(i),dp),aimag(zic(i))
            end do
            deallocate(bc,ac,yc,zic)
        case (3,4)
            allocate(w(nw))
            if (nw > 0) read(*,*) w
            if (operation == 3) then
                gd = group_delay(b,a,fs=fs,w=w)
                deallocate(b,a)
            else
                gd = group_delay(bc,ac,fs=fs,w=w)
                deallocate(bc,ac)
            end if
            print '(i0,1x,i0)',gd%status,size(gd%delay)
            if (size(gd%delay)>0) print '(*(es26.17e3,1x))',gd%delay
            deallocate(w)
        end select
    end do
end program filter_diagnostics_driver
