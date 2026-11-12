! SPDX-License-Identifier: MIT
! Stream independent complex SOS decimation outputs for SciPy comparison.
program reference_signal_decimate_complex_sos
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : decimate
    implicit none
    integer :: cases, index, mode, phase, n, q, ns, i, j, st
    real(dp), allocatable :: xr(:), xi(:), sr(:,:), si(:,:)
    complex(dp), allocatable :: x(:), sos(:,:), y(:)

    read(*,*) cases
    do index=1,cases
        read(*,*) mode,phase,n,q,ns
        allocate(xr(n),xi(n),sr(ns,6),si(ns,6),x(n),sos(ns,6))
        do j=1,ns
            read(*,*) sr(j,:)
            read(*,*) si(j,:)
        end do
        read(*,*) xr
        read(*,*) xi
        sos = cmplx(sr,si,dp)
        x = cmplx(xr,xi,dp)
        if (mode == 0) then
            y = decimate(xr,q,sos=sos,zero_phase=phase==1,status=st)
        else
            y = decimate(x,q,sos=sos,zero_phase=phase==1,status=st)
        end if
        print '(i0,1x,i0)',st,size(y)
        do i=1,size(y)
            print '(es26.17e3,1x,es26.17e3)',real(y(i),dp),aimag(y(i))
        end do
        deallocate(xr,xi,sr,si,x,sos,y)
    end do
end program reference_signal_decimate_complex_sos
