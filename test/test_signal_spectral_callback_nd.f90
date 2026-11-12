! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Callback-based detrending through assumed-rank spectral interfaces.
program test_signal_spectral_callback_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, &
        stft_result, spectrogram_result, welch, csd, coherence, periodogram, stft, spectrogram, &
        spectral_status_success, spectral_status_invalid_input
    use scifort_signal_spectral_complex, only : welch_complex, csd_complex, stft_complex, spectrogram_complex
    use scifort_signal_spectral_nd, only : welch_into, csd_into, coherence_into, periodogram_into
    use scifort_signal_timefreq_nd, only : stft_into, spectrogram_into
    implicit none
    real(dp), parameter :: tol = 5.0e-12_dp
    real(dp) :: x(64,3), y(64,3), w(10,3), cor(10,3), per(33,3)
    real(dp) :: strided_w(10,3), wrong(11,3), strided_output(20,3)
    real(dp) :: x15(32,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: out15(9,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: z(64,3), v(64,3), cross(10,3), cross_c(19,3)
    real(dp) :: x14(32,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: stft15(9,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    complex(dp) :: sx(10,3,7), sz(19,3,7), sp(10,3,7), spc(19,3,7)
    complex(dp) :: cs_sentinel(19,3,7)
    real(dp) :: spectral_bad(10,3), psd(10,3,7)
    real(dp), allocatable :: freq(:), times(:)
    type(spectral_result) :: r
    type(cross_spectral_result) :: c
    type(stft_result) :: f
    type(spectrogram_result) :: g
    integer :: i,j,stat, calls
    calls = 0
    do j=1,3
        do i=1,64
            x(i,j)=sin(0.13_dp*real(i*j,dp)) + 0.002_dp*real(i*i,dp) + 0.05_dp*j
            y(i,j)=cos(0.09_dp*real(i+j,dp)) + 0.003_dp*real(i*i,dp) - 0.11_dp*j
        end do
    end do
    z=cmplx(x,0.26_dp*y,dp)
    v=cmplx(y,0.19_dp*x,dp)

    w=-99.0_dp
    call welch_into(x,w,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=quadratic_real,status=stat)
    call valid(stat,'real welch into')
    if (size(freq)/=10) error stop 'wrong real welch frequency count'
    do j=1,3
        r=welch(x(:,j),nperseg=16,noverlap=8,nfft=19,detrend_callback=quadratic_real)
        call valid(r%status,'direct real welch')
        if (maxval(abs(w(:,j)-r%power))>tol) error stop 'real welch lane mismatch'
    end do
    call periodogram_into(x,per,freq,axis=1,detrend_callback=quadratic_real,status=stat)
    call valid(stat,'real periodogram into')
    do j=1,3
        r=periodogram(x(:,j),detrend_callback=quadratic_real)
        call valid(r%status,'direct periodogram')
        if (maxval(abs(per(:,j)-r%power))>tol) error stop 'periodogram lane mismatch'
    end do
    call coherence_into(x,y,cor,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=quadratic_real,status=stat)
    call valid(stat,'coherence into')
    do j=1,3
        r=coherence(x(:,j),y(:,j),nperseg=16,noverlap=8,nfft=19,detrend_callback=quadratic_real)
        call valid(r%status,'direct coherence')
        if (maxval(abs(cor(:,j)-r%power))>tol) error stop 'coherence lane mismatch'
    end do
    call csd_into(x,y,cross,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=quadratic_real,status=stat)
    call valid(stat,'real CSD into')
    do j=1,3
        c=csd(x(:,j),y(:,j),nperseg=16,noverlap=8,nfft=19,detrend_callback=quadratic_real)
        call valid(c%status,'direct real CSD')
        if (maxval(abs(cross(:,j)-c%power))>tol) error stop 'real CSD lane mismatch'
    end do
    call csd_into(z,v,cross_c,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=quadratic_complex,status=stat)
    call valid(stat,'complex CSD into')
    do j=1,3
        c=csd_complex(z(:,j),v(:,j),nperseg=16,noverlap=8,nfft=19,detrend_callback=quadratic_complex)
        call valid(c%status,'direct complex CSD')
        if (maxval(abs(cross_c(:,j)-c%power))>tol) error stop 'complex CSD lane mismatch'
    end do

    call csd_into(x,v,cross_c,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=quadratic_complex,status=stat)
    call valid(stat,'mixed real-complex CSD into')
    do j=1,3
        c=csd_complex(cmplx(x(:,j),0.0_dp,dp),v(:,j),nperseg=16,noverlap=8,nfft=19, &
            detrend_callback=quadratic_complex)
        call valid(c%status,'direct mixed CSD')
        if (maxval(abs(cross_c(:,j)-c%power))>tol) error stop 'mixed CSD lane mismatch'
    end do

    call stft_into(x,sx,freq,times,nperseg=16,noverlap=8,nfft=19, &
        boundary='none',padded=.false.,axis=1,detrend_callback=quadratic_real,status=stat)
    call valid(stat,'real stft into')
    if (size(freq)/=10 .or. size(times)/=7) error stop 'wrong STFT grids'
    do j=1,3
        f=stft(x(:,j),nperseg=16,noverlap=8,nfft=19,boundary='none',padded=.false., &
            detrend_callback=quadratic_real)
        call valid(f%status,'direct real stft')
        if (maxval(abs(sx(:,j,:)-f%values))>tol) error stop 'real STFT lane mismatch'
    end do
    call stft_into(z,sz,freq,times,nperseg=16,noverlap=8,nfft=19, &
        boundary='none',padded=.false.,axis=1,detrend_callback=quadratic_complex,status=stat)
    call valid(stat,'complex stft into')
    do j=1,3
        f=stft_complex(z(:,j),nperseg=16,noverlap=8,nfft=19,boundary='none',padded=.false., &
            detrend_callback=quadratic_complex)
        call valid(f%status,'direct complex stft')
        if (maxval(abs(sz(:,j,:)-f%values))>tol) error stop 'complex STFT lane mismatch'
    end do
    call spectrogram_into(x,psd,freq,times,nperseg=16,noverlap=8,nfft=19, &
        axis=1,mode='psd',detrend_callback=quadratic_real,status=stat)
    call valid(stat,'real PSD into')
    do j=1,3
        g=spectrogram(x(:,j),nperseg=16,noverlap=8,nfft=19, &
            mode='psd',detrend_callback=quadratic_real)
        call valid(g%status,'direct real PSD')
        if (maxval(abs(psd(:,j,:)-g%power))>tol) error stop 'real PSD lane mismatch'
    end do
    call spectrogram_into(x,sp,freq,times,nperseg=16,noverlap=8,nfft=19, &
        axis=1,mode='complex',detrend_callback=quadratic_real,status=stat)
    call valid(stat,'real complex-mode spectrogram into')
    do j=1,3
        g=spectrogram(x(:,j),nperseg=16,noverlap=8,nfft=19, &
            mode='complex',detrend_callback=quadratic_real)
        call valid(g%status,'direct complex-mode spectrogram')
        if (maxval(abs(sp(:,j,:)-g%values))>tol) error stop 'spectrogram lane mismatch'
    end do
    call spectrogram_into(z,spc,freq,times,nperseg=16,noverlap=8,nfft=19, &
        axis=1,mode='complex',detrend_callback=quadratic_complex,status=stat)
    call valid(stat,'complex spectrogram into')
    do j=1,3
        g=spectrogram_complex(z(:,j),nperseg=16,noverlap=8,nfft=19, &
            mode='complex',detrend_callback=quadratic_complex)
        call valid(g%status,'direct complex spectrogram')
        if (maxval(abs(spc(:,j,:)-g%values))>tol) error stop 'complex spectrogram lane mismatch'
    end do

    ! Signal arrays can be strided.  Verify no hidden contiguous assumption.
    call welch_into(x(1:64:2,:),strided_w,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=quadratic_real,status=stat)
    call valid(stat,'strided input')
    do j=1,3
        r=welch(x(1:64:2,j),nperseg=16,noverlap=8,nfft=19,detrend_callback=quadratic_real)
        call valid(r%status,'direct strided welch')
        if (maxval(abs(strided_w(:,j)-r%power))>tol) error stop 'strided lane mismatch'
    end do
    strided_output=-800.0_dp
    call welch_into(x,strided_output(1:20:2,:),freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=quadratic_real,status=stat)
    call valid(stat,'strided output')
    if (maxval(abs(strided_output(1:20:2,:)-w))>tol) error stop 'strided output values'
    if (any(strided_output(2:20:2,:)/=-800.0_dp)) error stop 'strided output surrounding values'
    wrong=111.0_dp
    call welch_into(x,wrong,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=quadratic_real,status=stat)
    if (stat==spectral_status_success .or. any(wrong/=111.0_dp)) error stop 'shape preservation failure'
    spectral_bad=123.0_dp
    call welch_into(x,spectral_bad,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend='constant',detrend_callback=quadratic_real,status=stat)
    if (stat==spectral_status_success .or. any(spectral_bad/=123.0_dp)) error stop 'ambiguous detrend rejection'
    spectral_bad=123.0_dp
    call welch_into(x,spectral_bad,freq,nperseg=16,noverlap=8,nfft=19,axis=1, &
        detrend_callback=nonfinite_real,status=stat)
    if (stat==spectral_status_success .or. any(spectral_bad/=123.0_dp)) error stop 'nonfinite real callback rejection'
    cs_sentinel=cmplx(7.0_dp,-3.0_dp,dp)
    call stft_into(z,cs_sentinel,freq,times,nperseg=16,noverlap=8,nfft=19, &
        boundary='none',padded=.false.,axis=1,detrend_callback=nonfinite_complex,status=stat)
    if (stat==spectral_status_success .or. any(abs(cs_sentinel-cmplx(7.0_dp,-3.0_dp,dp))>0.0_dp)) &
        error stop 'nonfinite complex callback rejection'

    ! Rank 15: one spectral lane and no intermediate rank-specific copy.
    x15=0.0_dp
    x15(:,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=x(1:32,1)
    out15=-1.0_dp
    call welch_into(x15,out15,freq,nperseg=16,noverlap=8,nfft=16,axis=1, &
        detrend_callback=quadratic_real,status=stat)
    call valid(stat,'rank 15 welch into')
    r=welch(x(1:32,1),nperseg=16,noverlap=8,nfft=16,detrend_callback=quadratic_real)
    if (maxval(abs(out15(:,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-r%power))>tol) &
        error stop 'rank15 spectral mismatch'
    x14=0.0_dp
    x14(:,1,1,1,1,1,1,1,1,1,1,1,1,1)=x(1:32,1)
    stft15=cmplx(-8.0_dp,-1.0_dp,dp)
    call stft_into(x14,stft15,freq,times,nperseg=16,noverlap=8,nfft=16, &
        boundary='none',padded=.false.,axis=1,detrend_callback=quadratic_real,status=stat)
    call valid(stat,'rank 14 STFT into rank 15')
    f=stft(x(1:32,1),nperseg=16,noverlap=8,nfft=16, &
        boundary='none',padded=.false.,detrend_callback=quadratic_real)
    if (maxval(abs(stft15(:,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-f%values))>tol) &
        error stop 'rank 14 STFT mismatch'
    print '(a,i0)', 'PASS native spectral callbacks, callback invocations: ', calls
contains
    subroutine valid(code,what)
        integer,intent(in) :: code !! zero is success
        character(*),intent(in) :: what !! failing context
        if (code/=spectral_status_success) then
            print '(a,a,1x,i0)', 'ERROR ',what,code
            error stop 'spectral callback error'
        end if
    end subroutine valid

    subroutine quadratic_real(samples)
        real(dp),intent(inout) :: samples(:) !! unwindowed samples
        real(dp) :: coordinate(size(samples)), sq(size(samples)), av, slope, curve
        integer :: k,n
        calls=calls+1
        n=size(samples)
        do k=1,n
            coordinate(k)=(real(k,dp)-0.5_dp*real(n+1,dp))/real(n,dp)
        end do
        sq=coordinate**2
        sq=sq-sum(sq)/real(n,dp)
        av=sum(samples)/real(n,dp)
        slope=sum(samples*coordinate)/sum(coordinate**2)
        curve=sum(samples*sq)/sum(sq**2)
        samples=samples-av-slope*coordinate-curve*sq
    end subroutine quadratic_real

    subroutine quadratic_complex(samples)
        complex(dp),intent(inout) :: samples(:) !! unwindowed complex samples
        real(dp) :: coordinate(size(samples)), sq(size(samples))
        complex(dp) :: av,slope,curve
        integer :: k,n
        calls=calls+1
        n=size(samples)
        do k=1,n
            coordinate(k)=(real(k,dp)-0.5_dp*real(n+1,dp))/real(n,dp)
        end do
        sq=coordinate**2
        sq=sq-sum(sq)/real(n,dp)
        av=sum(samples)/real(n,dp)
        slope=sum(samples*coordinate)/sum(coordinate**2)
        curve=sum(samples*sq)/sum(sq**2)
        samples=samples-av-slope*coordinate-curve*sq
    end subroutine quadratic_complex

    subroutine nonfinite_real(samples)
        real(dp),intent(inout) :: samples(:) !! deliberately invalid samples
        samples(1)=ieee_value(0.0_dp,ieee_quiet_nan)
    end subroutine nonfinite_real

    subroutine nonfinite_complex(samples)
        complex(dp),intent(inout) :: samples(:) !! deliberately invalid samples
        samples(1)=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    end subroutine nonfinite_complex
end program test_signal_spectral_callback_nd
