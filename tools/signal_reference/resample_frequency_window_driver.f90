! SPDX-License-Identifier: MIT
! Read batched independent signal.resample reference fixtures on standard input.
program resample_frequency_window_driver
    use scifort_kinds, only : dp
    use scifort_signal_multirate_nd_into, only : resample_into
    implicit none
    real(dp), allocatable :: xr(:,:), yr(:,:), wr(:)
    complex(dp), allocatable :: xc(:,:), yc(:,:), wc(:)
    integer :: n, m, nc, input_kind, window_kind, frequency_domain
    integer :: i, j, st, ios
    real(dp) :: ar, ai

    do
        read(*,*,iostat=ios) n,m,nc,input_kind,window_kind,frequency_domain
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid case header'
        allocate(xr(nc,n),xc(nc,n),yr(nc,m),yc(nc,m),wr(n),wc(n))
        do j=1,n
            read(*,*) ar, ai
            wr(j) = ar
            wc(j) = cmplx(ar,ai,dp)
        end do
        do j=1,n
            do i=1,nc
                read(*,*) ar, ai
                xr(i,j) = ar
                xc(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        if (input_kind == 0) then
            select case (window_kind)
            case (0)
                call resample_into(xr,yr,m,axis=2,frequency_window=wr,status=st)
            case (2)
                call resample_into(xr,yr,m,axis=2,window='hann',status=st)
            case default
                call resample_into(xr,yr,m,axis=2,status=st)
            end select
            yc = cmplx(yr,0.0_dp,dp)
        else if (window_kind == 2) then
            if (frequency_domain == 1) then
                call resample_into(xc,yc,m,axis=2,window='hann', &
                    domain='freq',status=st)
            else
                call resample_into(xc,yc,m,axis=2,window='hann',status=st)
            end if
        else if (window_kind == 3) then
            if (frequency_domain == 1) then
                call resample_into(xc,yc,m,axis=2,domain='freq',status=st)
            else
                call resample_into(xc,yc,m,axis=2,status=st)
            end if
        else if (window_kind == 0) then
            if (frequency_domain == 1) then
                call resample_into(xc,yc,m,axis=2,frequency_window=wr, &
                    domain='freq',status=st)
            else
                call resample_into(xc,yc,m,axis=2,frequency_window=wr,status=st)
            end if
        else
            if (frequency_domain == 1) then
                call resample_into(xc,yc,m,axis=2,complex_window=wc, &
                    domain='freq',status=st)
            else
                call resample_into(xc,yc,m,axis=2,complex_window=wc,status=st)
            end if
        end if
        print '(3(i0,1x))', st, nc, m
        do j=1,m
            do i=1,nc
                write(*,'(2(es26.17,1x))') real(yc(i,j),dp),aimag(yc(i,j))
            end do
        end do
        deallocate(xr,xc,yr,yc,wr,wc)
    end do
end program resample_frequency_window_driver
