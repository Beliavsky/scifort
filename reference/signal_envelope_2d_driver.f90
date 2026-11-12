! Reference-only standalone driver, not an FPM auto-test.
program signal_envelope_2d_driver
    use scifort_kinds, only : dp
    use scifort_signal_envelope, only : envelope
    implicit none
    integer :: nc, k, nr, nt, i, j, axis, nout, low, high, is_complex, square, residual, status
    real(dp) :: re, im
    real(dp), allocatable :: x(:,:), y(:,:), rr(:,:)
    complex(dp), allocatable :: z(:,:), cr(:,:)
    character(7) :: mode

    read(*,*) nc
    do k = 1, nc
        read(*,*) nr, nt, axis, nout, low, high, is_complex, square, residual
        allocate(x(nr,nt),z(nr,nt))
        do j = 1, nt
            do i = 1, nr
                read(*,*) re, im
                x(i,j) = re
                z(i,j) = cmplx(re,im,dp)
            end do
        end do
        select case(residual)
        case(0)
            mode = 'none'
        case(1)
            mode = 'all'
        case default
            mode = 'lowpass'
        end select
        if (is_complex == 0) then
            y = envelope(x, bp_low=low,bp_high=high,n_out=nout, squared=square==1, &
                residual=trim(mode),axis=axis,remainder=rr,status=status)
            write(*,'(5(i0,1x))') status, size(y,1), size(y,2), size(rr,1), size(rr,2)
            do j = 1, size(y,2)
                do i = 1, size(y,1)
                    re = 0.0_dp
                    if (size(rr,1)>0 .and. size(rr,2)>0) re = rr(i,j)
                    write(*,'(3(es26.17,1x))') y(i,j), re, 0.0_dp
                end do
            end do
        else
            y = envelope(z, bp_low=low,bp_high=high,n_out=nout,squared=square==1, &
                residual=trim(mode),axis=axis,remainder=cr,status=status)
            write(*,'(5(i0,1x))') status, size(y,1), size(y,2), size(cr,1), size(cr,2)
            do j = 1, size(y,2)
                do i = 1, size(y,1)
                    re = 0.0_dp
                    im = 0.0_dp
                    if (size(cr,1)>0 .and. size(cr,2)>0) then
                        re = real(cr(i,j),dp)
                        im = aimag(cr(i,j))
                    end if
                    write(*,'(3(es26.17,1x))') y(i,j), re, im
                end do
            end do
        end if
        deallocate(x,z,y)
        if (allocated(rr)) deallocate(rr)
        if (allocated(cr)) deallocate(cr)
    end do
end program signal_envelope_2d_driver
