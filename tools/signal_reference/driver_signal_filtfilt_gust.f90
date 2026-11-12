! SPDX-License-Identifier: MIT
! Reference interface for SciPy comparisons of Gustafsson filtering.
program driver_signal_filtfilt_gust
    use scifort_kinds, only : dp
    use scifort_signal_filtfilt, only : filtfilt
    use scifort_signal_filtfilt_gust, only : filtfilt_gust
    implicit none
    integer :: ncase, icase, mode, n, nb, na, irlen, i, status
    real(dp), allocatable :: br(:), ar(:), xr(:), yr(:)
    complex(dp), allocatable :: bc(:), ac(:), xc(:), yc(:)
    real(dp) :: re, im

    read(*,*) ncase
    do icase = 1, ncase
        read(*,*) mode, n, nb, na, irlen
        if (mode == 0) then
            allocate(br(nb), ar(na), xr(n))
            read(*,*) br
            read(*,*) ar
            read(*,*) xr
            if (irlen == -1) then
                yr = filtfilt(br, ar, xr, method='gust', status=status)
            else
                yr = filtfilt(br, ar, xr, method='gust', irlen=irlen, status=status)
            end if
            write(*,'(i0,1x,i0)') status, size(yr)
            do i = 1, size(yr)
                write(*,'(es26.17)') yr(i)
            end do
            deallocate(br,ar,xr,yr)
        else if (mode == 1) then
            allocate(br(nb), ar(na), xc(n))
            read(*,*) br
            read(*,*) ar
            do i=1,n
                read(*,*) re, im
                xc(i) = cmplx(re,im,dp)
            end do
            if (irlen == -1) then
                yc = filtfilt_gust(br, ar, xc, status=status)
            else
                yc = filtfilt_gust(br, ar, xc, irlen=irlen, status=status)
            end if
            write(*,'(i0,1x,i0)') status, size(yc)
            do i=1,size(yc)
                write(*,'(2(es26.17,1x))') real(yc(i),dp), aimag(yc(i))
            end do
            deallocate(br,ar,xc,yc)
        else
            allocate(bc(nb), ac(na), xc(n))
            do i=1,nb
                read(*,*) re, im
                bc(i) = cmplx(re,im,dp)
            end do
            do i=1,na
                read(*,*) re, im
                ac(i) = cmplx(re,im,dp)
            end do
            do i=1,n
                read(*,*) re, im
                xc(i) = cmplx(re,im,dp)
            end do
            if (irlen == -1) then
                yc = filtfilt(bc, ac, xc, method='gust', status=status)
            else
                yc = filtfilt(bc, ac, xc, method='gust', irlen=irlen, status=status)
            end if
            write(*,'(i0,1x,i0)') status, size(yc)
            do i=1,size(yc)
                write(*,'(2(es26.17,1x))') real(yc(i),dp), aimag(yc(i))
            end do
            deallocate(bc,ac,xc,yc)
        end if
    end do
end program driver_signal_filtfilt_gust
