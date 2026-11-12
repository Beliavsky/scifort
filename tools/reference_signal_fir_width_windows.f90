! SPDX-License-Identifier: MIT
! Deterministic SciFort FIR output cases for independent SciPy comparisons.
program reference_signal_fir_width_windows
    use scifort_kinds, only : dp
    use scifort_signal_design, only : firwin
    use scifort_signal_fir_advanced, only : firwin2
    implicit none
    real(dp), parameter :: fs_values(3) = [2.0_dp,40.0_dp,480.0_dp]
    real(dp), parameter :: width_ratios(4) = [0.004_dp,0.02_dp,0.10_dp,0.25_dp]
    real(dp), parameter :: gains_base(5) = [0.5_dp,0.90_dp,0.25_dp,0.60_dp,0.15_dp]
    real(dp) :: fs, nyq, width, gains(5), frequencies(5), cutoff(2)
    real(dp), allocatable :: taps(:), window_values(:)
    integer :: mode, ii, jj, kk, scale, n, nfreq, family, typ, st, grid
    logical :: anti
    character(8), parameter :: modes(4) = [character(8) :: &
        'lowpass', 'bandpass', 'highpass', 'bandstop']

    do mode=1,4
        do ii=1,2
            n = merge(15,26,ii==1)
            if (mode >= 3 .and. mod(n,2) == 0) n = n + 1
            do jj=1,size(fs_values)
                fs=fs_values(jj)
                nyq=fs/2.0_dp
                cutoff=[0.28_dp*nyq,0.62_dp*nyq]
                if (mode==2 .or. mode==4) cutoff(1)=0.25_dp*nyq
                do kk=1,size(width_ratios)
                    width=width_ratios(kk)*nyq
                    do scale=0,1
                        if (mode==1 .or. mode==3) then
                            taps=firwin(n,cutoff(:1),pass_zero=trim(modes(mode)), &
                                fs=fs,width=width,window='ignored-name', &
                                scale=scale==1,status=st)
                        else
                            taps=firwin(n,cutoff,pass_zero=trim(modes(mode)), &
                                fs=fs,width=width,window='ignored-name', &
                                scale=scale==1,status=st)
                        end if
                        call emit(taps, st)
                    end do
                end do
            end do
        end do
    end do

    do typ=1,4
        anti=typ==3 .or. typ==4
        do ii=1,2
            if (typ==1 .or. typ==3) then
                n=merge(13,29,ii==1)
            else
                n=merge(14,30,ii==1)
            end if
            fs=merge(2.0_dp,500.0_dp,ii==1)
            nyq=fs/2.0_dp
            frequencies=[0.0_dp,0.18_dp,0.5_dp,0.83_dp,1.0_dp]*nyq
            gains=gains_base
            if (anti) gains(1)=0.0_dp
            if (typ==2 .or. typ==3) gains(5)=0.0_dp
            do grid=1,2
                nfreq=merge(33,41,grid==1)
                if (n>=29) nfreq=merge(65,73,grid==1)
                do family=1,8
                    select case(family)
                    case(1)
                        taps=firwin2(n,frequencies,gains,nfreqs=nfreq, &
                            antisymmetric=anti,fs=fs,window='kaiser', &
                            window_parameter=5.2_dp,status=st)
                    case(2)
                        taps=firwin2(n,frequencies,gains,nfreqs=nfreq, &
                            antisymmetric=anti,fs=fs,window='gaussian', &
                            window_parameter=3.7_dp,status=st)
                    case(3)
                        taps=firwin2(n,frequencies,gains,nfreqs=nfreq, &
                            antisymmetric=anti,fs=fs,window='tukey', &
                            window_parameter=0.37_dp,status=st)
                    case(4)
                        taps=firwin2(n,frequencies,gains,nfreqs=nfreq, &
                            antisymmetric=anti,fs=fs,window='general_hamming', &
                            window_parameter=0.62_dp,status=st)
                    case(5)
                        taps=firwin2(n,frequencies,gains,nfreqs=nfreq, &
                            antisymmetric=anti,fs=fs,window='general_gaussian', &
                            window_parameter=1.6_dp,window_parameter2=4.4_dp,status=st)
                    case(6)
                        taps=firwin2(n,frequencies,gains,nfreqs=nfreq, &
                            antisymmetric=anti,fs=fs,window='general_cosine', &
                            window_coefficients=[0.4_dp,0.3_dp,0.2_dp,0.1_dp],status=st)
                    case(7)
                        taps=firwin2(n,frequencies,gains,nfreqs=nfreq, &
                            antisymmetric=anti,fs=fs,window='blackman',status=st)
                    case(8)
                        allocate(window_values(n))
                        window_values=[(0.2_dp+0.8_dp*sin(0.5_dp*acos(-1.0_dp)* &
                            real(kk-1,dp)/real(n-1,dp))**2,kk=1,n)]
                        taps=firwin2(n,frequencies,gains,nfreqs=nfreq, &
                            antisymmetric=anti,fs=fs,window_values=window_values,status=st)
                        deallocate(window_values)
                    end select
                    call emit(taps,st)
                end do
            end do
        end do
    end do
contains
    subroutine emit(weights, code)
        real(dp), intent(in) :: weights(:) !! output FIR coefficients
        integer, intent(in) :: code !! operation success status
        print '(i0,1x,i0)', code, size(weights)
        if (size(weights)>0) print '(*(es26.17,1x))', weights
    end subroutine emit
end program reference_signal_fir_width_windows
