! SPDX-License-Identifier: MIT
! Command-line protocol for reproducible analog IIR comparisons.
program analog_iir_driver
    use scifort_kinds, only : dp
    use scifort_signal_iir_extra, only : iirfilter
    use scifort_signal_design, only : butter_result
    implicit none
    type(butter_result) :: filt
    character(12) :: family, mode
    real(dp) :: attenuation, w1, w2
    integer :: order, ios, nb, na

    do
        read(*,*,iostat=ios) family, order, mode, attenuation, w1, w2
        if (ios < 0) exit
        if (ios /= 0) error stop 'analog IIR input parse'
        if (trim(mode) == 'bandpass' .or. trim(mode) == 'bandstop') then
            filt = iirfilter(order, [w1,w2], btype=trim(mode), ftype=trim(family), &
                rp=attenuation, rs=attenuation, analog=.true.)
        else
            filt = iirfilter(order, w1, btype=trim(mode), ftype=trim(family), &
                rp=attenuation, rs=attenuation, analog=.true.)
        end if
        nb = 0
        na = 0
        if (allocated(filt%b)) nb = size(filt%b)
        if (allocated(filt%a)) na = size(filt%a)
        write(*,'(i0,1x,i0,1x,i0)') filt%status,nb,na
        if (nb > 0) write(*,'(*(es25.16e3,1x))') filt%b
        if (na > 0) write(*,'(*(es25.16e3,1x))') filt%a
    end do
end program analog_iir_driver
