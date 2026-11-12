! SPDX-License-Identifier: MIT
! Batch test driver for analog and digital minimum-order IIR designs.
program order_band_driver
    use scifort_signal, only : dp, buttord, cheb1ord, cheb2ord, iirdesign, &
        iir_order_result, butter_result, design_status_success
    implicit none
    type(iir_order_result) :: sel
    type(butter_result) :: filt
    integer :: ios, family, amode, mode, nf, nb, na
    real(dp) :: gp, gs, fs, wp1, wp2, ws1, ws2, w1, w2
    logical :: analog
    character(7) :: ftype
    do
        read(*,*,iostat=ios) family, amode, mode, gp, gs, fs, wp1, wp2, ws1, ws2
        if (ios /= 0) exit
        analog = amode == 1
        select case (family)
        case (1)
            ftype = 'butter'
        case (2)
            ftype = 'cheby1'
        case (3)
            ftype = 'cheby2'
        case default
            error stop 'invalid family'
        end select
        if (mode <= 2) then
            if (analog) then
                select case (family)
                case (1)
                    sel = buttord(wp1,ws1,gp,gs,analog=.true.)
                case (2)
                    sel = cheb1ord(wp1,ws1,gp,gs,analog=.true.)
                case (3)
                    sel = cheb2ord(wp1,ws1,gp,gs,analog=.true.)
                end select
                filt = iirdesign(wp1,ws1,gp,gs,ftype=trim(ftype),analog=.true.)
            else
                select case (family)
                case (1)
                    sel = buttord(wp1,ws1,gp,gs,fs=fs)
                case (2)
                    sel = cheb1ord(wp1,ws1,gp,gs,fs=fs)
                case (3)
                    sel = cheb2ord(wp1,ws1,gp,gs,fs=fs)
                end select
                filt = iirdesign(wp1,ws1,gp,gs,ftype=trim(ftype),fs=fs)
            end if
        else
            if (analog) then
                select case (family)
                case (1)
                    sel = buttord([wp1,wp2],[ws1,ws2],gp,gs,analog=.true.)
                case (2)
                    sel = cheb1ord([wp1,wp2],[ws1,ws2],gp,gs,analog=.true.)
                case (3)
                    sel = cheb2ord([wp1,wp2],[ws1,ws2],gp,gs,analog=.true.)
                end select
                filt = iirdesign([wp1,wp2],[ws1,ws2],gp,gs,ftype=trim(ftype),analog=.true.)
            else
                select case (family)
                case (1)
                    sel = buttord([wp1,wp2],[ws1,ws2],gp,gs,fs=fs)
                case (2)
                    sel = cheb1ord([wp1,wp2],[ws1,ws2],gp,gs,fs=fs)
                case (3)
                    sel = cheb2ord([wp1,wp2],[ws1,ws2],gp,gs,fs=fs)
                end select
                filt = iirdesign([wp1,wp2],[ws1,ws2],gp,gs,ftype=trim(ftype),fs=fs)
            end if
        end if
        w1 = 0.0_dp
        w2 = 0.0_dp
        nf = 0
        if (allocated(sel%cutoffs)) then
            nf = size(sel%cutoffs)
            if (nf >= 1) w1 = sel%cutoffs(1)
            if (nf >= 2) w2 = sel%cutoffs(2)
        end if
        nb = 0
        na = 0
        if (allocated(filt%b)) nb = size(filt%b)
        if (allocated(filt%a)) na = size(filt%a)
        print '(3(i0,1x),2(es27.18,1x),3(i0,1x))', &
            sel%status, sel%order, nf, w1, w2, filt%status, nb, na
        if (filt%status == design_status_success) then
            print '(*(es27.18,1x))', filt%b
            print '(*(es27.18,1x))', filt%a
        end if
    end do
end program order_band_driver
