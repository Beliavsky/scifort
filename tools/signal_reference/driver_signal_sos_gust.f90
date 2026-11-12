! SPDX-License-Identifier: MIT
! Text protocol for independent Python SOS state-space Gustafsson comparisons.
program driver_signal_sos_gust
    use scifort_kinds, only : dp
    use scifort_signal_filtfilt_gust, only : sosfiltfilt_gust
    implicit none
    integer :: ncase, icase, mode, n, ns, irlen, i, j, code
    real(dp) :: re, im
    real(dp), allocatable :: sr(:,:), xr(:), yr(:)
    complex(dp), allocatable :: sc(:,:), xc(:), yc(:)

    read(*,*) ncase
    do icase = 1,ncase
        read(*,*) mode,n,ns,irlen
        allocate(sc(ns,6), xc(n))
        do i=1,ns
            do j=1,6
                read(*,*) re,im
                sc(i,j) = cmplx(re,im,dp)
            end do
        end do
        do i=1,n
            read(*,*) re,im
            xc(i) = cmplx(re,im,dp)
        end do
        select case(mode)
        case (0)
            sr = real(sc,dp)
            xr = real(xc,dp)
            if (irlen < 0) then
                yr = sosfiltfilt_gust(sr,xr,status=code)
            else
                yr = sosfiltfilt_gust(sr,xr,irlen=irlen,status=code)
            end if
            yc = cmplx(yr,0.0_dp,dp)
        case (1)
            sr = real(sc,dp)
            if (irlen < 0) then
                yc = sosfiltfilt_gust(sr,xc,status=code)
            else
                yc = sosfiltfilt_gust(sr,xc,irlen=irlen,status=code)
            end if
        case (2)
            xr = real(xc,dp)
            if (irlen < 0) then
                yc = sosfiltfilt_gust(sc,xr,status=code)
            else
                yc = sosfiltfilt_gust(sc,xr,irlen=irlen,status=code)
            end if
        case default
            if (irlen < 0) then
                yc = sosfiltfilt_gust(sc,xc,status=code)
            else
                yc = sosfiltfilt_gust(sc,xc,irlen=irlen,status=code)
            end if
        end select
        write(*,'(i0,1x,i0)') code,size(yc)
        do i=1,size(yc)
            write(*,'(2(es26.17,1x))') real(yc(i),dp),aimag(yc(i))
        end do
        deallocate(sc,xc)
    end do
end program driver_signal_sos_gust
