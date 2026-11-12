! SPDX-License-Identifier: MIT
! Numerical reference driver, input and output in two-line complex format.
program reference_signal_decimate_custom_filters
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : decimate
    implicit none
    integer :: cases, i, n, q, kind, mode, phase, nh, ns, j, status
    real(dp), allocatable :: realx(:), imagx(:), taps(:), sos(:,:), yr(:)
    complex(dp), allocatable :: z(:), yc(:)
    read(*,*) cases
    do i=1,cases
        read(*,*) kind,mode,phase,n,q,nh,ns
        allocate(realx(n),imagx(n),taps(nh),sos(ns,6),z(n))
        if (nh > 0) read(*,*) taps
        if (ns > 0) then
            do j=1,ns
                read(*,*) sos(j,:)
            end do
        end if
        read(*,*) realx
        read(*,*) imagx
        z = cmplx(realx,imagx,dp)
        if (kind == 0) then
            if (mode == 1) then
                yr = decimate(realx,q,taps=taps,zero_phase=phase==1,status=status)
            else
                yr = decimate(realx,q,sos=sos,zero_phase=phase==1,status=status)
            end if
            print '(i0,1x,i0)',status,size(yr)
            do j=1,size(yr)
                print '(es26.17e3,1x,es26.17e3)', yr(j),0.0_dp
            end do
            deallocate(yr)
        else
            if (mode == 1) then
                yc = decimate(z,q,taps=taps,zero_phase=phase==1,status=status)
            else
                yc = decimate(z,q,sos=sos,zero_phase=phase==1,status=status)
            end if
            print '(i0,1x,i0)',status,size(yc)
            do j=1,size(yc)
                print '(es26.17e3,1x,es26.17e3)',real(yc(j),dp),aimag(yc(j))
            end do
            deallocate(yc)
        end if
        deallocate(realx,imagx,taps,sos,z)
    end do
end program reference_signal_decimate_custom_filters
