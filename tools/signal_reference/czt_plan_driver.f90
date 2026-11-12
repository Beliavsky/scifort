! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Driver for deterministic comparisons against SciPy's CZT and ZoomFFT classes.
program czt_plan_driver
    use scifort_kinds, only : dp
    use scifort_signal_czt, only : czt_plan, zoom_fft_plan
    implicit none
    type(czt_plan) :: cz
    type(zoom_fft_plan) :: zz
    complex(dp), allocatable :: x(:,:), y(:,:)
    real(dp), allocatable :: xr(:), xi(:)
    real(dp) :: fn(2), fs, wr, wi, ar, ai
    integer :: cases, c, n, m, channels, axis, mode, real_input, endpoint, stat, nrow, ncol

    read(*,*) cases
    do c = 1,cases
        read(*,*) mode, n, m, channels, axis, real_input, endpoint
        read(*,*) fn, fs, wr, wi, ar, ai
        allocate(xr(n*channels),xi(n*channels))
        read(*,*) xr
        read(*,*) xi
        if (axis == 1) then
            allocate(x(n,channels))
        else
            allocate(x(channels,n))
        end if
        x = reshape(cmplx(xr,xi,dp),shape(x))
        if (mode <= 3) then
            select case(mode)
            case(1)
                call cz%initialize(n,m=m,status=stat)
            case default
                call cz%initialize(n,m=m,w=cmplx(wr,wi,dp),a=cmplx(ar,ai,dp),status=stat)
            end select
            if (stat == 0) then
                if (real_input == 1) then
                    y = cz%transform(real(x,dp),axis=axis,status=stat)
                else
                    y = cz%transform(x,axis=axis,status=stat)
                end if
            end if
        else
            if (mode == 4) then
                call zz%initialize(n,fn,m=m,fs=fs,endpoint=endpoint==1,status=stat)
            else
                call zz%initialize(n,fn(2),m=m,fs=fs,endpoint=endpoint==1,status=stat)
            end if
            if (stat == 0) then
                if (real_input == 1) then
                    y = zz%transform(real(x,dp),axis=axis,status=stat)
                else
                    y = zz%transform(x,axis=axis,status=stat)
                end if
            end if
        end if
        nrow = 0
        ncol = 0
        if (allocated(y)) then
            nrow = size(y,1)
            ncol = size(y,2)
        end if
        print '(3(i0,1x))', stat, nrow, ncol
        if (stat == 0) then
            print '(*(es25.16,1x))', real(y,dp)
            print '(*(es25.16,1x))', aimag(y)
        end if
        if (allocated(y)) deallocate(y)
        deallocate(x,xr,xi)
    end do
end program czt_plan_driver
