! SPDX-License-Identifier: MIT
! Independent SciPy-reference driver for parameterized and explicit windows.
program reference_signal_windowed_estimators
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : stft, spectrogram, welch, csd, coherence, &
        periodogram, stft_result, spectrogram_result, spectral_result, &
        cross_spectral_result, spectral_status_success
    use scifort_signal_spectral_complex, only : stft_complex, spectrogram_complex, &
        welch_complex, csd_complex, coherence_complex, periodogram_complex
    implicit none
    real(dp), parameter :: fs=7.25_dp
    real(dp) :: x(42), y(42)
    complex(dp) :: z(42), q(42)
    real(dp), allocatable :: p1,p2,coeff(:),weights(:),f(:),times(:),out(:)
    complex(dp), allocatable :: values(:)
    character(:), allocatable :: win
    integer :: method, style, dtype, i, state
    type(stft_result) :: t
    type(spectrogram_result) :: sp
    type(spectral_result) :: p
    type(cross_spectral_result) :: cp

    do i=1,42
        x(i) = sin(0.37_dp*real(i,dp)) + 0.012_dp*real(i,dp)
        y(i) = cos(0.11_dp*real(i,dp)) + 0.018_dp*real(i,dp) - 0.15_dp*x(i)
    end do
    z=cmplx(x,y,dp)
    q=cmplx(y,-0.25_dp*x,dp)
    do method=1,6
        do style=1,8
            do dtype=1,2
                if (allocated(win)) deallocate(win)
                if (allocated(p1)) deallocate(p1)
                if (allocated(p2)) deallocate(p2)
                if (allocated(coeff)) deallocate(coeff)
                if (allocated(weights)) deallocate(weights)
                select case(style)
                case(1)
                    win='tukey'
                    p1=0.32_dp
                case(2)
                    win='kaiser'
                    p1=6.2_dp
                case(3)
                    win='gaussian'
                    p1=2.8_dp
                case(4)
                    win='general_gaussian'
                    p1=1.4_dp
                    p2=3.0_dp
                case(5)
                    win='general_cosine'
                    coeff=[0.5_dp,0.4_dp,0.1_dp]
                case(6)
                    if (method == 6) then
                        allocate(weights(42))
                    else
                        allocate(weights(10))
                    end if
                    do i=1,size(weights)
                        weights(i)=0.6_dp+0.15_dp*sin(real(i,dp)*0.31_dp)
                    end do
                case(7)
                    ! Exercise estimator's default family (Tukey 0.25 in spectrogram).
                case(8)
                    win='tukey'
                    ! Plain named Tukey has alpha 0.5, unlike spectrogram default.
                end select
                allocate(f(0),times(0),values(0))
                select case(method)
                case(1) ! STFT complex spectrum / real one-sided spectrum
                    if (dtype==1) then
                        t=stft(x,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',boundary='odd',padded=.true.,scaling='psd', &
                            window_parameter=p1,window_parameter2=p2, &
                            window_coefficients=coeff,window_values=weights)
                    else
                        t=stft_complex(z,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',boundary='odd',padded=.true.,scaling='psd', &
                            window_parameter=p1,window_parameter2=p2, &
                            window_coefficients=coeff,window_values=weights)
                    end if
                    state=t%status
                    f=t%frequencies
                    times=t%times
                    values=reshape(t%values,[size(t%values)])
                case(2) ! spectrogram PSD
                    if (dtype==1) then
                        sp=spectrogram(x,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',scaling='density',window_parameter=p1, &
                            window_parameter2=p2,window_coefficients=coeff,window_values=weights)
                    else
                        sp=spectrogram_complex(z,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',scaling='density',window_parameter=p1, &
                            window_parameter2=p2,window_coefficients=coeff,window_values=weights)
                    end if
                    state=sp%status
                    f=sp%frequencies
                    times=sp%times
                    values=reshape(cmplx(sp%power,0.0_dp,dp),[size(sp%power)])
                case(3) ! Welch PSD
                    if (dtype==1) then
                        p=welch(x,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',scaling='density',window_parameter=p1, &
                            window_parameter2=p2,window_coefficients=coeff,window_values=weights)
                    else
                        p=welch_complex(z,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',scaling='density',window_parameter=p1, &
                            window_parameter2=p2,window_coefficients=coeff,window_values=weights)
                    end if
                    state=p%status
                    f=p%frequencies
                    values=cmplx(p%power,0.0_dp,dp)
                case(4) ! cross-spectral density
                    if (dtype==1) then
                        cp=csd(x,y,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',scaling='density',window_parameter=p1, &
                            window_parameter2=p2,window_coefficients=coeff,window_values=weights)
                    else
                        cp=csd_complex(z,q,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',scaling='density',window_parameter=p1, &
                            window_parameter2=p2,window_coefficients=coeff,window_values=weights)
                    end if
                    state=cp%status
                    f=cp%frequencies
                    values=cp%power
                case(5) ! coherence
                    if (dtype==1) then
                        p=coherence(x,y,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',window_parameter=p1,window_parameter2=p2, &
                            window_coefficients=coeff,window_values=weights)
                    else
                        p=coherence_complex(z,q,fs=fs,window=win,nperseg=10,noverlap=4,nfft=13, &
                            detrend='linear',window_parameter=p1,window_parameter2=p2, &
                            window_coefficients=coeff,window_values=weights)
                    end if
                    state=p%status
                    f=p%frequencies
                    values=cmplx(p%power,0.0_dp,dp)
                case(6) ! periodogram, 45 FFT points and full-length window
                    if (dtype==1) then
                        p=periodogram(x,fs=fs,window=win,nfft=45,detrend='linear', &
                            scaling='density',window_parameter=p1,window_parameter2=p2, &
                            window_coefficients=coeff,window_values=weights)
                    else
                        p=periodogram_complex(z,fs=fs,window=win,nfft=45,detrend='linear', &
                            scaling='density',window_parameter=p1,window_parameter2=p2, &
                            window_coefficients=coeff,window_values=weights)
                    end if
                    state=p%status
                    f=p%frequencies
                    values=cmplx(p%power,0.0_dp,dp)
                end select
                if (state/=spectral_status_success) then
                    print '(a,4(1x,i0))','FAILED',method,style,dtype,state
                    error stop 'window reference case failed'
                end if
                allocate(out(2*size(values)))
                out(1::2)=real(values,dp)
                out(2::2)=aimag(values)
                write(*,'(a,3(1x,i0))') 'CASE',method,style,dtype
                write(*,'(a,1x,i0,*(1x,es25.16))') 'FREQ',size(f),f
                write(*,'(a,1x,i0,*(1x,es25.16))') 'TIME',size(times),times
                write(*,'(a,1x,i0,*(1x,es25.16))') 'DATA',size(values),out
                deallocate(out,f,times,values)
            end do
        end do
    end do
end program reference_signal_windowed_estimators
