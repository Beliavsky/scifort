! SPDX-License-Identifier: MIT
! Numerical-reference producer for mixed real/complex CSD, coherence, and ShortTimeFFT.
program signal_mixed_cross_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, spectral_status_success
    use scifort_signal_spectral_complex, only : csd_real_complex, csd_complex_real, &
        coherence_real_complex, coherence_complex_real
    use scifort_signal_short_time_fft, only : short_time_fft, short_fft_success
    implicit none
    real(dp) :: x(79)
    complex(dp) :: y(79)
    complex(dp), allocatable :: values(:,:)
    type(cross_spectral_result) :: cross
    type(spectral_result) :: coh
    type(short_time_fft) :: plan
    character(24) :: window, trend, avg, scaling, mode, pad
    integer :: i, pair, wid, length_id, overlap_id, nfft_id, avg_id, trend_id, scale_id
    integer :: nw, no, nf, k, p, mode_id, pad_id, phase_id, status, hop, phase

    do i=1,size(x)
        x(i) = 0.7_dp*cos(0.19_dp*i) + sin(0.43_dp*i) + 0.013_dp*i
        y(i) = cmplx(0.4_dp*sin(0.23_dp*i) + x(i), &
            cos(0.37_dp*i) - 0.15_dp*sin(0.09_dp*i), dp)
    end do
    do pair=1,2
        do wid=1,2
            window = 'hann'
            if (wid == 2) window = 'gaussian'
            do length_id=1,2
                nw = 12 + 3*(length_id-1)
                do overlap_id=1,2
                    no = 3 + 4*(overlap_id-1)
                    do nfft_id=1,2
                        nf = 19 + 5*(nfft_id-1)
                        do avg_id=1,2
                            avg = 'mean'
                            if (avg_id == 2) avg = 'median'
                            do trend_id=1,2
                                trend = 'none'
                                if (trend_id == 2) trend = 'linear'
                                do scale_id=1,2
                                    scaling = 'density'
                                    if (scale_id == 2) scaling = 'spectrum'
                                    if (pair == 1) then
                                        if (wid == 1) then
                                            cross = csd_real_complex(x,y,fs=20.0_dp,window=trim(window), &
                                                nperseg=nw,noverlap=no,nfft=nf,detrend=trim(trend), &
                                                scaling=trim(scaling),average=trim(avg))
                                        else
                                            cross = csd_real_complex(x,y,fs=20.0_dp,window=trim(window), &
                                                window_parameter=3.5_dp,nperseg=nw,noverlap=no,nfft=nf, &
                                                detrend=trim(trend),scaling=trim(scaling),average=trim(avg))
                                        end if
                                    else
                                        if (wid == 1) then
                                            cross = csd_complex_real(y,x,fs=20.0_dp,window=trim(window), &
                                                nperseg=nw,noverlap=no,nfft=nf,detrend=trim(trend), &
                                                scaling=trim(scaling),average=trim(avg))
                                        else
                                            cross = csd_complex_real(y,x,fs=20.0_dp,window=trim(window), &
                                                window_parameter=3.5_dp,nperseg=nw,noverlap=no,nfft=nf, &
                                                detrend=trim(trend),scaling=trim(scaling),average=trim(avg))
                                        end if
                                    end if
                                    if (cross%status /= spectral_status_success) error stop 'mixed CSD reference'
                                    write(*,'(a,1x,10(i0,1x))') 'C', pair, wid, nw, no, nf, avg_id, trend_id, scale_id, &
                                        size(cross%power), cross%status
                                    do k=1,size(cross%power)
                                        write(*,'(3(es26.17,1x))') cross%frequencies(k), &
                                            real(cross%power(k),dp), aimag(cross%power(k))
                                    end do
                                end do
                            end do
                        end do
                        do trend_id=1,2
                            trend = 'none'
                            if (trend_id == 2) trend = 'linear'
                            if (pair == 1) then
                                if (wid == 1) then
                                    coh = coherence_real_complex(x,y,fs=20.0_dp,window=trim(window), &
                                        nperseg=nw,noverlap=no,nfft=nf,detrend=trim(trend))
                                else
                                    coh = coherence_real_complex(x,y,fs=20.0_dp,window=trim(window), &
                                        window_parameter=3.5_dp,nperseg=nw,noverlap=no,nfft=nf, &
                                        detrend=trim(trend))
                                end if
                            else
                                if (wid == 1) then
                                    coh = coherence_complex_real(y,x,fs=20.0_dp,window=trim(window), &
                                        nperseg=nw,noverlap=no,nfft=nf,detrend=trim(trend))
                                else
                                    coh = coherence_complex_real(y,x,fs=20.0_dp,window=trim(window), &
                                        window_parameter=3.5_dp,nperseg=nw,noverlap=no,nfft=nf, &
                                        detrend=trim(trend))
                                end if
                            end if
                            if (coh%status /= spectral_status_success) error stop 'mixed coherence reference'
                            write(*,'(a,1x,8(i0,1x))') 'H',pair,wid,nw,no,nf,trend_id,size(coh%power),coh%status
                            do k=1,size(coh%power)
                                write(*,'(2(es26.17,1x))') coh%frequencies(k),coh%power(k)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    do pair=1,2
        do mode_id=1,2
            mode = 'twosided'
            if (mode_id == 2) mode = 'centered'
            do pad_id=1,4
                select case(pad_id)
                case(1)
                    pad = 'zeros'
                case(2)
                    pad = 'edge'
                case(3)
                    pad = 'even'
                case default
                    pad = 'odd'
                end select
                do phase_id=1,2
                    phase = 2*(phase_id-1)
                    do scale_id=1,2
                        scaling = 'psd'
                        if (scale_id == 2) scaling = 'magnitude'
                        nw = 9
                        nf = 13
                        hop = 3
                        call plan%initialize_window('hann',nw,hop,20.0_dp, &
                            mfft=nf,fft_mode=trim(mode),scale_to=trim(scaling),phase_shift=phase)
                        if (plan%status /= short_fft_success) error stop 'mixed ShortTimeFFT init'
                        if (pair == 1) then
                            call plan%cross_spectrogram(x,y,values,status,padding=trim(pad))
                        else
                            call plan%cross_spectrogram(y,x,values,status,padding=trim(pad))
                        end if
                        if (status /= short_fft_success) error stop 'mixed ShortTimeFFT result'
                        write(*,'(a,1x,11(i0,1x))') 'S',pair,mode_id,pad_id,phase,scale_id, &
                            nw,hop,nf,size(x),size(values,1),size(values,2)
                        do p=1,size(values,2)
                            do k=1,size(values,1)
                                write(*,'(2(es26.17,1x))') real(values(k,p),dp),aimag(values(k,p))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end program signal_mixed_cross_reference_driver
