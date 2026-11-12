! SPDX-License-Identifier: MIT
! Standalone stream driver for comparisons with SciPy's custom dlti decimate.
program reference_signal_decimate_transfer_function
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : decimate
    implicit none
    integer :: cases, index, mode, phase, n, q, nb, na, i, st
    real(dp), allocatable :: xr(:), xi(:), br(:), bi(:), ar(:), ai(:), out_real(:)
    complex(dp), allocatable :: x(:), b(:), a(:), y(:)
    read(*,*) cases
    do index = 1, cases
        read(*,*) mode,phase,n,q,nb,na
        allocate(xr(n),xi(n),br(nb),bi(nb),ar(na),ai(na),x(n),b(nb),a(na))
        read(*,*) br
        read(*,*) bi
        read(*,*) ar
        read(*,*) ai
        read(*,*) xr
        read(*,*) xi
        x = cmplx(xr,xi,dp)
        b = cmplx(br,bi,dp)
        a = cmplx(ar,ai,dp)
        select case(mode)
        case (0)
            out_real = decimate(xr,q,b=br,a=ar,zero_phase=phase==1,status=st)
            print '(i0,1x,i0)', st,size(out_real)
            do i=1,size(out_real)
                print '(es26.17e3,1x,es26.17e3)',out_real(i),0.0_dp
            end do
        case (1)
            y = decimate(x,q,b=b,a=a,zero_phase=phase==1,status=st)
            print '(i0,1x,i0)',st,size(y)
            do i=1,size(y)
                print '(es26.17e3,1x,es26.17e3)',real(y(i),dp),aimag(y(i))
            end do
        case (2)
            y = decimate(xr,q,b=b,a=a,zero_phase=phase==1,status=st)
            print '(i0,1x,i0)',st,size(y)
            do i=1,size(y)
                print '(es26.17e3,1x,es26.17e3)',real(y(i),dp),aimag(y(i))
            end do
        case default
            error stop 'bad comparison mode'
        end select
        deallocate(xr,xi,br,bi,ar,ai,x,b,a)
        if (allocated(y)) deallocate(y)
        if (allocated(out_real)) deallocate(out_real)
    end do
end program reference_signal_decimate_transfer_function
