! SPDX-License-Identifier: MIT
! Fortran numerical-reference producer for compare_signal_short_time_fft.py
program test_signal_short_time_fft_reference
    use scifort_kinds, only : dp
    use scifort_signal_short_time_fft, only : short_time_fft, short_time_fft_result, &
        short_time_spectrogram_result, short_fft_success
    implicit none
    type(short_time_fft) :: plan
    type(short_time_fft_result) :: spectrum
    type(short_time_spectrogram_result) :: powers
    real(dp), allocatable :: x(:), xr(:)
    complex(dp), allocatable :: z(:), zr(:), cross(:,:)
    character(24) :: mode, scale, pad, kind, trend
    integer :: i, ic, n, nw, hop, nf, phase, q, p, ncases, status, p0, p1, offset
    real(dp) :: fs

    open(unit=11,file='short_fft_cases.txt',status='old',action='read')
    open(unit=12,file='short_fft_results.txt',status='replace',action='write')
    read(11,*) ncases
    do ic=1,ncases
        read(11,*) n,nw,hop,nf,fs,mode,scale,phase,pad,kind,p0,p1,offset,trend
        allocate(x(n),z(n))
        do i=1,n
            x(i) = cos(0.17_dp*real(i-1,dp))+0.3_dp*sin(0.32_dp*real(i-1,dp))
            z(i) = cmplx(x(i),0.37_dp*cos(0.09_dp*real(i-1,dp)),dp)
        end do
        call plan%initialize_window('hann',nw,hop,fs,mfft=nf,fft_mode=trim(mode), &
            scale_to=trim(scale),phase_shift=phase)
        if (plan%status /= short_fft_success) error stop 'bad reference init'
        if (trim(kind) == 'real') then
            if (trim(trend) == 'none') then
                call plan%stft(x,spectrum,p0=p0,p1=p1,k_offset=offset,padding=trim(pad))
            else
                call plan%stft_detrend(x,spectrum,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad))
            end if
        else
            if (trim(trend) == 'none') then
                call plan%stft(z,spectrum,p0=p0,p1=p1,k_offset=offset,padding=trim(pad))
            else
                call plan%stft_detrend(z,spectrum,trim(trend),p0=p0,p1=p1, &
                    k_offset=offset,padding=trim(pad))
            end if
        end if
        if (spectrum%status /= short_fft_success) error stop 'bad reference stft'
        call plan%cross_spectrogram(x,x,cross,status,padding=trim(pad))
        if (status /= short_fft_success .and. trim(kind)=='real') error stop 'cross status'
        if (p0 == plan%p_min .and. p1 == plan%p_max(n) .and. &
            offset == 0 .and. trim(trend) == 'none') then
        if (trim(mode)=='onesided' .or. trim(mode)=='onesided2X') then
            call plan%istft(spectrum%values,xr,status,k1=n)
            if (status /= short_fft_success) error stop 'bad reference real inverse'
            if (maxval(abs(xr-x)) > 2.0e-10_dp) then
                print *, 'fail',ic,n,nw,hop,nf,trim(mode),trim(scale),phase,trim(pad),maxval(abs(xr-x))
                error stop 'bad reference real reconstruction'
            end if
        else
            call plan%istft(spectrum%values,zr,status,k1=n)
            if (status /= short_fft_success) error stop 'bad reference complex inverse'
            if (trim(kind)=='real') then
                if (maxval(abs(zr-cmplx(x,0.0_dp,dp))) > 2.0e-10_dp) then
                    print *, 'fail complex from real',ic,n,nw,hop,nf,trim(mode),trim(scale),phase,trim(pad),maxval(abs(zr-cmplx(x,0.0_dp,dp)))
                    error stop 'bad real-to-complex reconstruction'
                end if
            else
                if (maxval(abs(zr-z)) > 2.0e-10_dp) then
                    print *, 'fail complex',ic,n,nw,hop,nf,trim(mode),trim(scale),phase,trim(pad),maxval(abs(zr-z))
                    error stop 'bad reference complex reconstruction'
                end if
            end if
        end if
        end if
        write(12,*) ic, size(spectrum%values,1),size(spectrum%values,2),spectrum%p0
        do p=1,size(spectrum%values,2)
            do q=1,size(spectrum%values,1)
                write(12,'(2(es26.17,1x))') real(spectrum%values(q,p),dp),aimag(spectrum%values(q,p))
            end do
        end do
        deallocate(x,z)
    end do
    close(12)
    close(11)
    print *, 'short-time FFT reference cases written:', ncases
end program test_signal_short_time_fft_reference
