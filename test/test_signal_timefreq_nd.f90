! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_timefreq_nd
    use scifort_signal, only : dp, stft_into, spectrogram_into, &
        stft, spectrogram, stft_result, spectrogram_result, &
        spectral_status_success, spectral_status_invalid_input
    implicit none
    real(dp) :: x(8,2), x2(2,8), x3(2,8,2), xr14(1,1,1,1,1,1,1,1,1,1,1,1,2,8)
    real(dp) :: p(5,2,3), p3(2,5,2,3), &
        p15(1,1,1,1,1,1,1,1,1,1,1,1,2,5,3)
    real(dp) :: strided(8,4), ps(5,4,3), bad(4,2,3), phases(5,2,3)
    complex(dp) :: z(8,2), zs(8,2,3), sp(5,2,3), sp3(2,5,2,3)
    complex(dp) :: spec(8,2,3), sp15(1,1,1,1,1,1,1,1,1,1,1,1,2,5,3)
    real(dp), allocatable :: f(:),t(:)
    type(stft_result) :: s
    type(spectrogram_result) :: q
    integer :: i,j,k,st
    real(dp), parameter :: tol=2.0e-11_dp

    do j=1,2
        do i=1,8
            x(i,j)=sin(0.19_dp*i+0.1_dp*j)+0.2_dp*cos(0.31_dp*i-0.4_dp*j)
        end do
    end do
    z=cmplx(x,0.15_dp*x,dp)
    sp=(-300.0_dp,20.0_dp)
    call stft_into(x,sp,f,t,nperseg=4,noverlap=2,nfft=8, &
        boundary='none',padded=.false.,axis=1,status=st)
    if(st/=spectral_status_success .or. size(f)/=5 .or. size(t)/=3) &
        error stop 'stft rank2 real status'
    do j=1,2
        s=stft(x(:,j),nperseg=4,noverlap=2,nfft=8,boundary='none',padded=.false.)
        if(maxval(abs(s%values-sp(:,j,:)))>tol) error stop 'stft rank2 real mismatch'
    end do

    p=-1.0_dp
    call spectrogram_into(x,p,f,t,nperseg=4,noverlap=2,nfft=8,axis=1,status=st)
    if(st/=spectral_status_success) error stop 'spectrogram status'
    do j=1,2
        q=spectrogram(x(:,j),nperseg=4,noverlap=2,nfft=8)
        if(maxval(abs(q%power-p(:,j,:)))>tol) error stop 'spectrogram rank2 mismatch'
    end do
    spec=(-300.0_dp,20.0_dp)
    call spectrogram_into(z,spec,f,t,nperseg=4,noverlap=2,nfft=8,axis=-2, &
        mode='complex',status=st)
    if(st/=spectral_status_success .or. size(f)/=8) error stop 'spectrogram complex status'
    do j=1,2
        s=stft(z(:,j),window='tukey',nperseg=4,noverlap=2,nfft=8, &
            detrend='constant',boundary='none',padded=.false.,scaling='psd')
        if(maxval(abs(s%values-spec(:,j,:)))>tol) error stop 'spectrogram complex mismatch'
    end do
    phases=-200.0_dp
    call spectrogram_into(x,phases,f,t,nperseg=4,noverlap=2,nfft=8, &
        axis=1,mode='phase',status=st)
    if(st/=spectral_status_success .or. any(phases==-200.0_dp)) error stop 'phase result'
    bad=14.0_dp
    call spectrogram_into(x,bad,f,t,nperseg=4,noverlap=2,nfft=8,axis=1,status=st)
    if(st/=spectral_status_invalid_input .or. any(bad/=14.0_dp)) error stop 'bad shape unchanged'
    sp=(-30.0_dp,20.0_dp)
    call stft_into(x,sp,f,t,nperseg=4,noverlap=4,nfft=8,axis=1,status=st)
    if(st/=spectral_status_invalid_input .or. any(sp/=cmplx(-30.0_dp,20.0_dp,dp))) &
        error stop 'bad option unchanged'
    call stft_into(x,sp,f,t,nperseg=4,noverlap=2,nfft=8,axis=0,status=st)
    if(st/=spectral_status_invalid_input .or. any(sp/=cmplx(-30.0_dp,20.0_dp,dp))) &
        error stop 'bad axis unchanged'
    call spectrogram_into(x,sp,f,t,nperseg=4,noverlap=2,nfft=8,axis=1, &
        mode='psd',status=st)
    if(st/=spectral_status_invalid_input .or. any(sp/=cmplx(-30.0_dp,20.0_dp,dp))) &
        error stop 'wrong output type unchanged'
    do j=1,2
        x2(j,:)=x(:,j)
        do k=1,2
            x3(k,:,j)=x(:,j)+0.02_dp*k
        end do
    end do
    sp3=(-50.0_dp,10.0_dp)
    call stft_into(x3,sp3,f,t,nperseg=4,noverlap=2,nfft=8, &
        boundary='none',padded=.false.,axis=2,status=st)
    if(st/=spectral_status_success) error stop 'rank3 stft status'
    do j=1,2
        do k=1,2
            s=stft(x3(k,:,j),nperseg=4,noverlap=2,nfft=8,boundary='none',padded=.false.)
            if(maxval(abs(s%values-sp3(k,:,j,:)))>tol) error stop 'rank3 mismatch'
        end do
    end do
    strided(:,1:4:2)=x
    strided(:,2:4:2)=0.0_dp
    ps=-100.0_dp
    call spectrogram_into(strided(:,1:4:2),ps(:,2:4:2,:),f,t, &
        nperseg=4,noverlap=2,nfft=8,axis=1,status=st)
    if(st/=spectral_status_success) error stop 'strided status'
    if(any(ps(:,1:4:2,:)/=-100.0_dp)) error stop 'strided unrelated destination'
    do j=1,2
        q=spectrogram(x(:,j),nperseg=4,noverlap=2,nfft=8)
        if(maxval(abs(q%power-ps(:,2*j,:)))>tol) error stop 'strided output mismatch'
    end do
    xr14=0.0_dp
    xr14(1,1,1,1,1,1,1,1,1,1,1,1,1,:)=x(:,1)
    xr14(1,1,1,1,1,1,1,1,1,1,1,1,2,:)=x(:,2)
    p15=0.0_dp
    call spectrogram_into(xr14,p15,f,t,nperseg=4,noverlap=2,nfft=8,axis=-1,status=st)
    if(st/=spectral_status_success) error stop 'rank14 spectrogram status'
    do j=1,2
        q=spectrogram(x(:,j),nperseg=4,noverlap=2,nfft=8)
        if(maxval(abs(q%power-p15(1,1,1,1,1,1,1,1,1,1,1,1,j,:,:)))>tol) &
            error stop 'rank14 spectrogram mismatch'
    end do
    sp15=cmplx(-20.0_dp,0.0_dp,dp)
    call stft_into(xr14,sp15,f,t,nperseg=4,noverlap=2,nfft=8, &
        boundary='none',padded=.false.,axis=-1,status=st)
    if(st/=spectral_status_success) error stop 'rank14 stft status'
    if(maxval(abs(sp15(1,1,1,1,1,1,1,1,1,1,1,1,1,:,:)- &
        stft_values(x(:,1))))>tol) error stop 'rank14 stft mismatch'
    print '(a)', 'test_signal_timefreq_nd passed'
contains
    function stft_values(v) result(arr)
        real(dp), intent(in) :: v(:) !! test input
        complex(dp), allocatable :: arr(:,:)
        type(stft_result) :: s0
        s0=stft(v,nperseg=4,noverlap=2,nfft=8,boundary='none',padded=.false.)
        arr=s0%values
    end function stft_values
end program test_signal_timefreq_nd
