! SPDX-License-Identifier: MIT
! Test-only batch driver for SciPy elliptic IIR comparisons.
program elliptic_driver
    use scifort_kinds, only : dp
    use scifort_signal_iir_extra, only : ellip, ellipord, iirdesign, iir_order_result
    use scifort_signal_design, only : butter_result
    implicit none
    integer :: ios, order, mode, analog_flag, edges
    real(dp) :: rp, rs, fs, cutoff(2), wp(2), ws(2)
    type(butter_result) :: filt
    type(iir_order_result) :: selection
    character(len=1) :: op
    character(len=9) :: btype

    do
        read(*,*,iostat=ios) op, order, mode, analog_flag, rp, rs, fs, cutoff, wp, ws
        if (ios /= 0) exit
        select case (mode)
        case (1)
            btype = 'lowpass'
        case (2)
            btype = 'highpass'
        case (3)
            btype = 'bandpass'
        case (4)
            btype = 'bandstop'
        case default
            error stop 'invalid test mode'
        end select
        edges = 1
        if (mode > 2) edges = 2
        select case (op)
        case ('F')
            if (analog_flag /= 0) then
                filt = ellip(order,rp,rs,cutoff(:edges),btype=btype,analog=.true.)
            else
                filt = ellip(order,rp,rs,cutoff(:edges),btype=btype,fs=fs)
            end if
        case ('D')
            if (analog_flag /= 0) then
                filt = iirdesign(wp(:edges),ws(:edges),rp,rs,ftype='ellip',analog=.true.)
            else
                filt = iirdesign(wp(:edges),ws(:edges),rp,rs,ftype='ellip',fs=fs)
            end if
        case ('O')
            if (analog_flag /= 0) then
                selection = ellipord(wp(:edges),ws(:edges),rp,rs,analog=.true.)
            else
                selection = ellipord(wp(:edges),ws(:edges),rp,rs,fs=fs)
            end if
            print '(i0,1x,i0,1x,i0)',selection%status,selection%order, &
                allocated_size(selection%cutoffs)
            if (selection%status == 0) print '(*(es27.18,1x))',selection%cutoffs
            cycle
        case default
            error stop 'unknown elliptic comparison operation'
        end select
        print '(i0,1x,i0,1x,i0)',filt%status,allocated_size(filt%b), &
            allocated_size(filt%a)
        if (filt%status == 0) then
            print '(*(es27.18,1x))',filt%b
            print '(*(es27.18,1x))',filt%a
        end if
    end do
contains
    integer function allocated_size(values) result(n)
        real(dp), allocatable, intent(in) :: values(:) !! optional allocated array
        n = 0
        if (allocated(values)) n = size(values)
    end function allocated_size
end program elliptic_driver
