! SPDX-License-Identifier: MIT
! Deterministic ShortTimeFFT callable detrending output for SciPy comparison.
program signal_stft_callable_detrend_driver
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result, short_fft_success
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result) :: spectrum
    real(dp) :: x(48)
    complex(dp) :: z(48)
    character(12) :: mode, padding, scaling, window
    integer :: i, typ, wid, imode, ipad, iscale, iphase
    integer :: nw, hop, nf, phase, j, k

    do i=1,size(x)
        x(i) = sin(0.15_dp*i)+0.2_dp*cos(0.4_dp*i)+0.003_dp*i*i
        z(i) = cmplx(x(i),0.1_dp*cos(0.3_dp*i),dp)
    end do
    do typ=1,2
        do wid=1,2
            nw=7+wid
            hop=2+wid
            nf=nw+5
            window='hann'
            if (wid==2) window='hamming'
            do imode=1,4
                if (typ==2 .and. imode<3) cycle
                select case (imode)
                case(1)
                    mode='onesided'
                case(2)
                    mode='onesided2X'
                case(3)
                    mode='twosided'
                case(4)
                    mode='centered'
                end select
                do ipad=1,4
                    select case (ipad)
                    case(1)
                        padding='zeros'
                    case(2)
                        padding='edge'
                    case(3)
                        padding='even'
                    case(4)
                        padding='odd'
                    end select
                    do iscale=1,3
                        scaling='none'
                        if (iscale==2) scaling='magnitude'
                        if (iscale==3) scaling='psd'
                        if (imode==2 .and. iscale==1) cycle
                        do iphase=1,2
                            phase=2*(iphase-1)
                            call plan%initialize_window(trim(window),nw,hop,20.0_dp, &
                                mfft=nf,fft_mode=trim(mode),scale_to=trim(scaling),phase_shift=phase)
                            if (plan%status /= short_fft_success) error stop 'plan init'
                            if (typ==1) then
                                call plan%stft_detrend(x,spectrum,detrend=nonlinear_real,padding=trim(padding))
                            else
                                call plan%stft_detrend(z,spectrum,detrend=nonlinear_complex,padding=trim(padding))
                            end if
                            if (spectrum%status /= short_fft_success) error stop 'callback STFT'
                            write(*,'(a,1x,11(i0,1x))') 'S',typ,wid,nw,hop,nf,imode,ipad,iscale,phase, &
                                size(spectrum%values,1),size(spectrum%values,2)
                            write(*,'(i0)') spectrum%p0
                            do j=1,size(spectrum%frequencies)
                                write(*,'(es26.17)') spectrum%frequencies(j)
                            end do
                            do j=1,size(spectrum%times)
                                write(*,'(es26.17)') spectrum%times(j)
                            end do
                            do j=1,size(spectrum%values,2)
                                do k=1,size(spectrum%values,1)
                                    write(*,'(2(es26.17,1x))') real(spectrum%values(k,j),dp), &
                                        aimag(spectrum%values(k,j))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
contains
    subroutine nonlinear_real(samples)
        real(dp), intent(inout) :: samples(:) !! one real STFT segment
        samples = samples - sum(samples)/real(size(samples),dp) + 0.07_dp*samples**2
    end subroutine nonlinear_real
    subroutine nonlinear_complex(samples)
        complex(dp), intent(inout) :: samples(:) !! one complex STFT segment
        samples = samples - sum(samples)/real(size(samples),dp) + 0.07_dp*samples**2 &
            + cmplx(0.0_dp,0.04_dp,dp)*samples
    end subroutine nonlinear_complex
end program signal_stft_callable_detrend_driver
