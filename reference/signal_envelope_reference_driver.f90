! SPDX-License-Identifier: MIT
! Reference driver compiled separately, not an FPM auto-test.
program signal_envelope_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_envelope, only : envelope
    implicit none
    integer :: cases, i, j, n, nout, lo, hi, cx, squared, method, status
    complex(dp), allocatable :: z(:), rcomplex(:)
    real(dp), allocatable :: xr(:), yr(:), rr(:)
    real(dp) :: realpart, imagpart
    character(7) :: name

    read(*,*) cases
    do i = 1, cases
        read(*,*) n, nout, lo, hi, cx, squared, method
        allocate(z(n),xr(n))
        do j = 1, n
            read(*,*) realpart, imagpart
            z(j) = cmplx(realpart,imagpart,dp)
            xr(j) = realpart
        end do
        select case (method)
        case (0)
            name = 'none'
        case (1)
            name = 'all'
        case default
            name = 'lowpass'
        end select
        if (cx == 0) then
            yr = envelope(xr, bp_low=lo, bp_high=hi, n_out=nout, squared=(squared==1), &
                residual=trim(name), remainder=rr, status=status)
            write(*,'(3(i0,1x))') status, size(yr), size(rr)
            do j = 1, size(yr)
                if (j <= size(rr)) then
                    write(*,'(3(es25.16,1x))') yr(j), rr(j), 0.0_dp
                else
                    write(*,'(3(es25.16,1x))') yr(j), 0.0_dp, 0.0_dp
                end if
            end do
        else
            yr = envelope(z, bp_low=lo, bp_high=hi, n_out=nout, squared=(squared==1), &
                residual=trim(name), remainder=rcomplex, status=status)
            write(*,'(3(i0,1x))') status, size(yr), size(rcomplex)
            do j = 1, size(yr)
                if (j <= size(rcomplex)) then
                    write(*,'(3(es25.16,1x))') yr(j), real(rcomplex(j),dp), aimag(rcomplex(j))
                else
                    write(*,'(3(es25.16,1x))') yr(j), 0.0_dp, 0.0_dp
                end if
            end do
        end if
        deallocate(z,xr)
        if (allocated(rr)) deallocate(rr)
        if (allocated(rcomplex)) deallocate(rcomplex)
        if (allocated(yr)) deallocate(yr)
    end do
end program signal_envelope_reference_driver
