! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! stdin: numtaps ncut pass_zero window_id scale_flag width_flag scalar_flag fs
!        followed by 5 cutoff positions (only the first ncut are used).
! stdout: status count [count coefficients], one record per input line.
program reference_signal_firwin_multiband
    use, intrinsic :: iso_fortran_env, only : input_unit
    use scifort_kinds, only : dp
    use scifort_signal_design, only : firwin
    implicit none
    integer :: numtaps, ncut, pass_id, window_id, scale_flag, width_flag, scalar_flag, &
        ios, st, i
    real(dp) :: fs, bounds(5)
    real(dp), allocatable :: taps(:)
    character(:), allocatable :: window, passmode
    logical :: scale_val, logical_pass

    do
        read(input_unit,*,iostat=ios) numtaps,ncut,pass_id,window_id, &
            scale_flag,width_flag,scalar_flag,fs,bounds
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid reference input record'
        if (ncut < 1 .or. ncut > 5) error stop 'invalid reference cutoff count'
        select case (window_id)
        case (0)
            window='hamming'
        case (1)
            window='hann'
        case (2)
            window='blackman'
        case (3)
            window='kaiser'
        case default
            error stop 'invalid window id'
        end select
        select case (pass_id)
        case (2)
            passmode='lowpass'
        case (3)
            passmode='highpass'
        case (4)
            passmode='bandpass'
        case (5)
            passmode='bandstop'
        case default
            passmode=''
        end select
        scale_val = scale_flag /= 0
        logical_pass = pass_id == 1
        ! Test scalar/string, scalar/logical, vector/string, vector/logical,
        ! and omitted pass_zero resolution with the same randomized records.
        if (scalar_flag == 1) then
            select case (pass_id)
            case (0,1)
                if (window_id == 3) then
                    if (width_flag == 1) then
                        taps = firwin(numtaps,bounds(1),pass_zero=logical_pass, &
                            fs=fs,scale=scale_val,width=0.04_dp*fs,status=st)
                    else
                        taps = firwin(numtaps,bounds(1),pass_zero=logical_pass, &
                            fs=fs,scale=scale_val,window=window, &
                            window_parameter=3.2_dp,status=st)
                    end if
                else
                    if (width_flag == 1) then
                        taps = firwin(numtaps,bounds(1),pass_zero=logical_pass, &
                            fs=fs,scale=scale_val,width=0.04_dp*fs,status=st)
                    else
                        taps = firwin(numtaps,bounds(1),pass_zero=logical_pass, &
                            fs=fs,scale=scale_val,window=window,status=st)
                    end if
                end if
            case (6)
                if (width_flag == 1) then
                    taps = firwin(numtaps,bounds(1),fs=fs,scale=scale_val, &
                        width=0.04_dp*fs,status=st)
                else if (window_id == 3) then
                    taps = firwin(numtaps,bounds(1),fs=fs,scale=scale_val, &
                        window=window,window_parameter=3.2_dp,status=st)
                else
                    taps = firwin(numtaps,bounds(1),fs=fs,scale=scale_val, &
                        window=window,status=st)
                end if
            case default
                if (width_flag == 1) then
                    taps = firwin(numtaps,bounds(1),pass_zero=passmode, &
                        fs=fs,scale=scale_val,width=0.04_dp*fs,status=st)
                else if (window_id == 3) then
                    taps = firwin(numtaps,bounds(1),pass_zero=passmode, &
                        fs=fs,scale=scale_val,window=window, &
                        window_parameter=3.2_dp,status=st)
                else
                    taps = firwin(numtaps,bounds(1),pass_zero=passmode, &
                        fs=fs,scale=scale_val,window=window,status=st)
                end if
            end select
        else
            select case (pass_id)
            case (0,1)
                if (width_flag == 1) then
                    taps = firwin(numtaps,bounds(:ncut),pass_zero=logical_pass, &
                        fs=fs,scale=scale_val,width=0.04_dp*fs,status=st)
                else if (window_id == 3) then
                    taps = firwin(numtaps,bounds(:ncut),pass_zero=logical_pass, &
                        fs=fs,scale=scale_val,window=window, &
                        window_parameter=3.2_dp,status=st)
                else
                    taps = firwin(numtaps,bounds(:ncut),pass_zero=logical_pass, &
                        fs=fs,scale=scale_val,window=window,status=st)
                end if
            case (6)
                if (width_flag == 1) then
                    taps = firwin(numtaps,bounds(:ncut),fs=fs,scale=scale_val, &
                        width=0.04_dp*fs,status=st)
                else if (window_id == 3) then
                    taps = firwin(numtaps,bounds(:ncut),fs=fs,scale=scale_val, &
                        window=window,window_parameter=3.2_dp,status=st)
                else
                    taps = firwin(numtaps,bounds(:ncut),fs=fs,scale=scale_val, &
                        window=window,status=st)
                end if
            case default
                if (width_flag == 1) then
                    taps = firwin(numtaps,bounds(:ncut),pass_zero=passmode, &
                        fs=fs,scale=scale_val,width=0.04_dp*fs,status=st)
                else if (window_id == 3) then
                    taps = firwin(numtaps,bounds(:ncut),pass_zero=passmode, &
                        fs=fs,scale=scale_val,window=window, &
                        window_parameter=3.2_dp,status=st)
                else
                    taps = firwin(numtaps,bounds(:ncut),pass_zero=passmode, &
                        fs=fs,scale=scale_val,window=window,status=st)
                end if
            end select
        end if
        write(*,'(i0,1x,i0,*(1x,es24.16))') st,size(taps),(taps(i),i=1,size(taps))
    end do
end program reference_signal_firwin_multiband
