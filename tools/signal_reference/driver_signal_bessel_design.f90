! SPDX-License-Identifier: MIT
! Numerical reference driver for SciFort Bessel IIR design.
program driver_signal_bessel_design
    use scifort_kinds, only : dp
    use scifort_signal_bessel_design, only : bessel
    use scifort_signal_design, only : butter_result, design_status_success
    implicit none
    type(butter_result) :: filt
    integer :: ncase, caseidx, order, mode_id, analog_id, norm_id, ios, nsec
    real(dp) :: fs, edge1, edge2
    character(8), parameter :: modes(4) = [character(8) :: &
        'lowpass','highpass','bandpass','bandstop']
    character(5), parameter :: norms(3) = [character(5) :: &
        'phase','delay','mag']

    read(*,*,iostat=ios) ncase
    if (ios /= 0) error stop 'case count'
    do caseidx = 1, ncase
        read(*,*,iostat=ios) order,mode_id,analog_id,norm_id,fs,edge1,edge2
        if (ios /= 0) error stop 'case parameters'
        if (mode_id <= 2) then
            if (analog_id == 1) then
                filt = bessel(order,edge1,btype=trim(modes(mode_id)), &
                    norm=trim(norms(norm_id)),analog=.true.)
            else
                filt = bessel(order,edge1,btype=trim(modes(mode_id)), &
                    norm=trim(norms(norm_id)),fs=fs)
            end if
        else
            if (analog_id == 1) then
                filt = bessel(order,[edge1,edge2],btype=trim(modes(mode_id)), &
                    norm=trim(norms(norm_id)),analog=.true.)
            else
                filt = bessel(order,[edge1,edge2],btype=trim(modes(mode_id)), &
                    norm=trim(norms(norm_id)),fs=fs)
            end if
        end if
        if (filt%status /= design_status_success) then
            print '(i0,1x,i0,1x,i0,1x,i0)', filt%status,0,0,0
            cycle
        end if
        nsec = 0
        if (allocated(filt%sos)) nsec = size(filt%sos,1)
        print '(i0,1x,i0,1x,i0,1x,i0)', filt%status, size(filt%b), size(filt%a), nsec
        print '(*(es25.16,1x))', filt%b
        print '(*(es25.16,1x))', filt%a
        if (nsec > 0) then
            ! Use section-major flattening for Python reshaping.
            print '(*(es25.16,1x))', transpose(filt%sos)
        end if
    end do
end program driver_signal_bessel_design
