! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Input/output driver for compare_signal_analog_sos_scipy.py (not an fpm test).
program analog_sos_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_analog_sos, only : analog_zpk_result, analog_sos_result, &
        iirfilter_zpk, iirfilter_sos, sosfreqs
    use scifort_signal_design, only : frequency_response
    implicit none
    character(8), parameter :: shapes(4) = ['lowpass ', 'highpass', 'bandpass', 'bandstop']
    character(7), parameter :: families(4) = ['butter ', 'cheby1 ', 'cheby2 ', 'ellip  ']
    type(analog_zpk_result) :: zpk
    type(analog_sos_result) :: sections
    type(frequency_response) :: response
    real(dp) :: lo, hi, rp, rs, w(10)
    real(dp), allocatable :: row(:)
    integer :: order, family, shape, ios, i, idx, nz, np

    do
        read(*,*,iostat=ios) order, family, shape, lo, hi, rp, rs
        if (ios /= 0) exit
        if (family < 1 .or. family > 4 .or. shape < 1 .or. shape > 4) stop 1
        if (shape <= 2) then
            zpk = iirfilter_zpk(order,lo,btype=trim(shapes(shape)), &
                ftype=trim(families(family)),rp=rp,rs=rs)
            sections = iirfilter_sos(order,lo,btype=trim(shapes(shape)), &
                ftype=trim(families(family)),rp=rp,rs=rs)
        else
            zpk = iirfilter_zpk(order,[lo,hi],btype=trim(shapes(shape)), &
                ftype=trim(families(family)),rp=rp,rs=rs)
            sections = iirfilter_sos(order,[lo,hi],btype=trim(shapes(shape)), &
                ftype=trim(families(family)),rp=rp,rs=rs)
        end if
        if (zpk%status /= 0 .or. sections%status /= 0) then
            write(*,'(a)') 'FAILED'
            cycle
        end if
        w = [0.0_dp,lo*0.1_dp,lo*0.3_dp,lo*0.8_dp,lo,lo*1.2_dp, &
            hi*0.8_dp,hi,hi*1.2_dp,hi*10.0_dp]
        response = sosfreqs(sections%sos,w)
        if (response%status /= 0) then
            write(*,'(a)') 'FAILED'
            cycle
        end if
        nz = size(zpk%zeros)
        np = size(zpk%poles)
        allocate(row(2*(nz+np+size(w))))
        idx = 0
        do i = 1, nz
            row(idx+1:idx+2) = [real(zpk%zeros(i),dp), aimag(zpk%zeros(i))]
            idx = idx+2
        end do
        do i = 1, np
            row(idx+1:idx+2) = [real(zpk%poles(i),dp), aimag(zpk%poles(i))]
            idx = idx+2
        end do
        do i = 1, size(w)
            row(idx+1:idx+2) = [real(response%response(i),dp),aimag(response%response(i))]
            idx = idx+2
        end do
        write(*,'(i0,1x,i0,1x,es25.17)') nz,np,zpk%gain
        write(*,'(*(es25.17,1x))') row
        deallocate(row)
    end do
end program analog_sos_reference_driver
