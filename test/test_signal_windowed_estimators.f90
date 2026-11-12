! SPDX-License-Identifier: MIT
! Verify window parameter/vector forwarding, dimensions, and transactional output.
program test_signal_windowed_estimators
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : stft, spectrogram, welch, csd, &
        coherence, periodogram, stft_result, spectrogram_result, &
        spectral_result, cross_spectral_result, spectral_status_success, &
        spectral_status_invalid_input
    use scifort_signal_spectral_complex, only : stft_complex, spectrogram_complex, &
        welch_complex, csd_complex, periodogram_complex
    use scifort_signal_spectral_nd, only : periodogram_into, welch_into, &
        csd_into, coherence_into
    use scifort_signal_timefreq_nd, only : stft_into, spectrogram_into
    implicit none
    real(dp), parameter :: tol=3.0e-11_dp
    real(dp) :: x(12,2),y(12,2),weights(6),long_window(12),coeff(3)
    real(dp) :: pxx(8,2), pgram(7,2), coh(5,2), psd(5,2,3), magnitude(8,2,3)
    real(dp) :: xr14(1,1,1,1,1,1,1,1,1,1,1,1,1,12)
    real(dp) :: pr15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,7)
    complex(dp) :: zx(12,2), zy(12,2), cxy(5,2), sreal(5,2,5), &
        scomplex(8,2,5), spec_complex(5,2,3), sparse(5,4,5)
    complex(dp) :: sc15(1,1,1,1,1,1,1,1,1,1,1,1,1,5,5)
    real(dp) :: array_noncontig(12,4)
    real(dp), allocatable :: frequencies(:),times(:)
    type(stft_result) :: st, stc
    type(spectral_result) :: sr
    type(cross_spectral_result) :: cr
    type(spectrogram_result) :: spec
    integer :: i,j,code

    weights = [(0.65_dp+0.025_dp*real(i,dp),i=1,6)]
    long_window = [(0.5_dp+0.015_dp*real(i,dp),i=1,12)]
    coeff = [0.5_dp,0.4_dp,0.1_dp]
    do j=1,2
        do i=1,12
            x(i,j)=sin(0.24_dp*real(i*j,dp))+0.023_dp*real(i,dp)
            y(i,j)=cos(0.2_dp*real(i*j,dp))+0.1_dp*x(i,j)
        end do
    end do
    zx=cmplx(x,0.34_dp*y,dp)
    zy=cmplx(y,-0.13_dp*x,dp)

    ! A real window vector fixes nperseg even if the caller omits that keyword.
    st=stft(x(:,1),window_values=weights,noverlap=3,nfft=8,boundary='odd', &
        padded=.false.,detrend='none')
    if (st%status/=spectral_status_success) error stop 'explicit STFT window inference'
    if (size(st%values,1)/=5 .or. size(st%values,2)/=5) &
        error stop 'explicit STFT lengths'
    sreal=(99.0_dp,-11.0_dp)
    call stft_into(x,sreal,frequencies,times,window_values=weights,noverlap=3, &
        nfft=8,boundary='odd',padded=.false.,detrend='none',axis=1,status=code)
    if (code/=spectral_status_success) error stop 'real multidimensional STFT'
    do j=1,2
        st=stft(x(:,j),window_values=weights,noverlap=3,nfft=8, &
            boundary='odd',padded=.false.,detrend='none')
        if (maxval(abs(sreal(:,j,:)-st%values))>tol) error stop 'real STFT values'
    end do

    scomplex=(99.0_dp,11.0_dp)
    call stft_into(zx,scomplex,frequencies,times,window='gaussian', &
        window_parameter=2.8_dp,nperseg=6,noverlap=3,nfft=8, &
        boundary='even',padded=.false.,detrend='linear',axis=-2,status=code)
    if (code/=spectral_status_success) error stop 'complex multidimensional STFT'
    do j=1,2
        stc=stft_complex(zx(:,j),window='gaussian',window_parameter=2.8_dp, &
            nperseg=6,noverlap=3,nfft=8,boundary='even',padded=.false.,detrend='linear')
        if (maxval(abs(scomplex(:,j,:)-stc%values))>tol) &
            error stop 'complex STFT values'
    end do

    call spectrogram_into(x,psd,frequencies,times,nperseg=6,noverlap=3,nfft=8, &
        detrend='constant',axis=1,status=code)
    if (code/=spectral_status_success) error stop 'default Tukey spectrogram'
    do j=1,2
        spec=spectrogram(x(:,j),nperseg=6,noverlap=3,nfft=8)
        if (maxval(abs(psd(:,j,:)-spec%power))>tol) error stop 'spectrogram PSD'
    end do
    call spectrogram_into(x,spec_complex,frequencies,times,nperseg=6, &
        noverlap=3,nfft=8,mode='complex',axis=1,status=code)
    if (code/=spectral_status_success) error stop 'complex spectrogram mode'
    do j=1,2
        st=stft(x(:,j),window='tukey',window_parameter=0.25_dp,nperseg=6, &
            noverlap=3,nfft=8,detrend='constant',boundary='none', &
            padded=.false.,scaling='psd')
        if (maxval(abs(spec_complex(:,j,:)-st%values))>tol) &
            error stop 'complex spectrogram default alpha'
    end do
    call spectrogram_into(zx,magnitude,frequencies,times,window_values=weights, &
        noverlap=3,nfft=8,mode='magnitude',axis=1,status=code)
    if (code/=spectral_status_success) error stop 'magnitude explicit window'
    do j=1,2
        stc=stft_complex(zx(:,j),window_values=weights,noverlap=3,nfft=8, &
            detrend='constant',boundary='none',padded=.false.,scaling='psd')
        if (maxval(abs(magnitude(:,j,:)-abs(stc%values)))>tol) &
            error stop 'complex magnitude values'
    end do

    ! Periodogram / Welch / CSD / coherence assumed-rank paths.
    call periodogram_into(x,pgram,frequencies,window_values=long_window, &
        nfft=12,axis=1,status=code)
    if (code/=spectral_status_success) error stop 'periodogram caller window'
    do j=1,2
        sr=periodogram(x(:,j),window_values=long_window,nfft=12)
        if (maxval(abs(pgram(:,j)-sr%power))>tol) error stop 'periodogram values'
    end do
    call welch_into(zx,pxx,frequencies,window='kaiser',window_parameter=6.0_dp, &
        nperseg=6,noverlap=3,nfft=8,axis=1,status=code)
    if (code/=spectral_status_success) error stop 'Welch tuple window'
    do j=1,2
        sr=welch_complex(zx(:,j),window='kaiser',window_parameter=6.0_dp, &
            nperseg=6,noverlap=3,nfft=8)
        if (maxval(abs(pxx(:,j)-sr%power))>tol) &
            error stop 'Welch values'
    end do
    call csd_into(x,y,cxy,frequencies,window='general_cosine', &
        window_coefficients=coeff,nperseg=6,noverlap=3,nfft=8, &
        axis=1,status=code)
    if (code/=spectral_status_success) error stop 'CSD cosine window'
    do j=1,2
        cr=csd(x(:,j),y(:,j),window='general_cosine', &
            window_coefficients=coeff,nperseg=6,noverlap=3,nfft=8)
        if (maxval(abs(cxy(:,j)-cr%power))>tol) error stop 'CSD values'
    end do
    call coherence_into(x,y,coh,frequencies,window_values=weights, &
        noverlap=3,nfft=8,axis=1,status=code)
    if (code/=spectral_status_success) error stop 'coherence caller window'
    do j=1,2
        sr=coherence(x(:,j),y(:,j),window_values=weights,noverlap=3,nfft=8)
        if (maxval(abs(coh(:,j)-sr%power))>tol) error stop 'coherence values'
    end do

    ! Error checking remains atomic: bad arguments do not modify output arrays.
    sparse=(123.0_dp,-56.0_dp)
    call stft_into(x,sparse(:,1:4:2,:),frequencies,times,window='hann', &
        window_values=weights,nperseg=6,noverlap=3,nfft=8,axis=1,status=code)
    if (code/=spectral_status_invalid_input) error stop 'incompatible windows accepted'
    if (any(sparse/=(123.0_dp,-56.0_dp))) error stop 'error modified complex output'
    if (size(frequencies)/=0 .or. size(times)/=0) error stop 'failure coordinate vectors'
    call stft_into(x,sparse(:,1:4:2,:),frequencies,times,window_values=weights, &
        nperseg=5,noverlap=3,nfft=8,axis=1,status=code)
    if (code/=spectral_status_invalid_input) error stop 'window length mismatch accepted'
    if (any(sparse/=(123.0_dp,-56.0_dp))) error stop 'error modified complex output'
    st=stft(x(:,1),window='hann',window_parameter=0.45_dp, &
        nperseg=6,noverlap=3,nfft=8)
    if (st%status/=spectral_status_invalid_input) error stop 'unused window parameter accepted'
    st=stft(x(:,1),window='kaiser',window_parameter=3.0_dp, &
        window_parameter2=4.0_dp,nperseg=6)
    if (st%status/=spectral_status_invalid_input) error stop 'unused second parameter accepted'
    long_window(7)=ieee_value(0.0_dp,ieee_quiet_nan)
    pgram=-9.0_dp
    call periodogram_into(x,pgram,frequencies,window_values=long_window, &
        nfft=12,axis=1,status=code)
    if (code/=spectral_status_invalid_input) error stop 'nonfinite window accepted'
    if (any(pgram/=-9.0_dp)) error stop 'error modified periodogram'
    long_window(7)=1.0_dp

    ! Strided input and output sections must be handled without aliasing.
    array_noncontig=-800.0_dp
    array_noncontig(:,1:4:2)=x
    sparse=(123.0_dp,-56.0_dp)
    call stft_into(array_noncontig(:,1:4:2),sparse(:,1:4:2,:),frequencies,times, &
        window_values=weights,noverlap=3,nfft=8,detrend='none',boundary='odd', &
        padded=.false.,axis=1,status=code)
    if (code/=spectral_status_success) error stop 'strided STFT'
    if (maxval(abs(sparse(:,1:4:2,:)-sreal))>tol) error stop 'strided STFT values'
    if (any(sparse(:,2:4:2,:)/=(123.0_dp,-56.0_dp))) &
        error stop 'strided STFT affected unused cells'

    ! One sample dimension in rank 14 becomes frequency dimension in rank 15.
    xr14(1,1,1,1,1,1,1,1,1,1,1,1,1,:)=x(:,1)
    sc15=(-10.0_dp,12.0_dp)
    call stft_into(xr14,sc15,frequencies,times,window_values=weights, &
        noverlap=3,nfft=8,boundary='odd',padded=.false.,detrend='none', &
        axis=-1,status=code)
    if (code/=spectral_status_success) error stop 'rank14 explicit-window STFT'
    st=stft(x(:,1),window_values=weights,noverlap=3,nfft=8, &
        boundary='odd',padded=.false.,detrend='none')
    if (maxval(abs(sc15(1,1,1,1,1,1,1,1,1,1,1,1,1,:,:)-st%values))>tol) &
        error stop 'rank14 STFT values'
    pr15=123.0_dp
    block
        real(dp) :: xr15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,12)
        xr15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)=x(:,1)
        call periodogram_into(xr15,pr15,frequencies,window_values=long_window, &
            nfft=12,axis=-1,status=code)
        if (code/=spectral_status_success) error stop 'rank15 caller-window periodogram'
    end block
    sr=periodogram(x(:,1),window_values=long_window,nfft=12)
    if (maxval(abs(pr15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-sr%power))>tol) &
        error stop 'rank15 periodogram values'

    print '(a)','test_signal_windowed_estimators passed'
end program test_signal_windowed_estimators
