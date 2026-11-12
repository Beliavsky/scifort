! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_spectral_nd
    use scifort_signal, only : dp, periodogram, welch, csd, coherence, &
        periodogram_into, welch_into, csd_into, coherence_into, spectral_status_success, &
        spectral_status_invalid_input, spectral_result, cross_spectral_result
    implicit none
    real(dp) :: x(8,3), y(8,3), expected(5,3), power(5,3), p3(2,5,3)
    real(dp) :: x3(2,8,3), y3(2,8,3), short(4,3), wrong(3,3)
    real(dp) :: padded_y(4,3), strided(8,6), strided_output(5,6)
    real(dp) :: xr15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,8)
    real(dp) :: ps15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,5)
    complex(dp) :: z(8,3), cz(8,3), cross(8,3), xmix(2,8,3), cross3(2,8,3)
    complex(dp) :: c15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,8)
    complex(dp) :: cp15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,8)
    real(dp), allocatable :: f(:)
    type(spectral_result) :: p
    type(cross_spectral_result) :: c
    integer :: i, j, k, st
    real(dp), parameter :: tol = 3.0e-11_dp

    do j=1,3
        do i=1,8
            x(i,j) = sin(real(i+j,dp)*0.35_dp) + 0.1_dp*real(j,dp)
            y(i,j) = cos(real(2*i+j,dp)*0.18_dp)
        end do
    end do
    z = cmplx(x,0.1_dp*y,dp)
    cz = cmplx(y,-0.2_dp*x,dp)
    power = -123.0_dp
    call periodogram_into(x,power,f,axis=1,fs=4.0_dp,status=st)
    if(st/=spectral_status_success .or. size(f)/=5) error stop 'periodogram into status'
    do j=1,3
        p=periodogram(x(:,j),fs=4.0_dp)
        if(maxval(abs(p%power-power(:,j)))>tol) error stop 'periodogram lane'
    end do
    power = -1.0_dp
    call welch_into(x,power,f,axis=1,fs=4.0_dp,window='hann', &
        nperseg=4,noverlap=2,nfft=8,average='median',status=st)
    if(st/=spectral_status_success) error stop 'welch median status'
    do j=1,3
        p=welch(x(:,j),fs=4.0_dp,nperseg=4,noverlap=2,nfft=8,average='median')
        if(maxval(abs(p%power-power(:,j)))>tol) error stop 'welch median lane'
    end do
    cross = (-77.0_dp,-33.0_dp)
    call csd_into(z,cz,cross,f,axis=1,nperseg=4,noverlap=2,nfft=8,status=st)
    if(st/=spectral_status_success .or. size(f)/=8) error stop 'complex csd status'
    do j=1,3
        c=csd(z(:,j),cz(:,j),nperseg=4,noverlap=2,nfft=8)
        if(maxval(abs(c%power-cross(:,j)))>tol) error stop 'complex csd lane'
    end do
    cross = (-77.0_dp,-33.0_dp)
    call csd_into(x,cz,cross,f,axis=1,nperseg=4,noverlap=2,nfft=8,status=st)
    if(st/=spectral_status_success) error stop 'mixed csd status'
    do j=1,3
        c=csd(cmplx(x(:,j),0.0_dp,dp),cz(:,j),nperseg=4,noverlap=2,nfft=8)
        if(maxval(abs(c%power-cross(:,j)))>tol) error stop 'mixed csd lane'
    end do
    power=-5.0_dp
    call coherence_into(x,y,power,f,axis=1,nperseg=4,noverlap=2,nfft=8,status=st)
    if(st/=spectral_status_success) error stop 'coherence into status'
    do j=1,3
        p=coherence(x(:,j),y(:,j),nperseg=4,noverlap=2,nfft=8)
        if(maxval(abs(p%power-power(:,j)))>tol) error stop 'coherence lane'
    end do
    do j=1,3
        do k=1,2
            x3(k,:,j)=x(:,j)+0.03_dp*k
            y3(k,:,j)=y(:,j)-0.1_dp*k
        end do
    end do
    p3=-100.0_dp
    call periodogram_into(x3,p3,f,axis=2,status=st)
    if(st/=spectral_status_success) error stop 'rank3 status'
    do j=1,3
        do k=1,2
            p=periodogram(x3(k,:,j))
            if(maxval(abs(p%power-p3(k,:,j)))>tol) error stop 'rank3 lane'
        end do
    end do
    xmix=cmplx(x3,0.2_dp*y3,dp)
    cross3=(-100.0_dp,12.0_dp)
    call csd_into(x3,xmix,cross3,f,axis=2,nperseg=4,noverlap=2,nfft=8,status=st)
    if(st/=spectral_status_success) error stop 'rank3 mixed status'
    do j=1,3
        do k=1,2
            c=csd(cmplx(x3(k,:,j),0.0_dp,dp),xmix(k,:,j),nperseg=4,noverlap=2,nfft=8)
            if(maxval(abs(c%power-cross3(k,:,j)))>tol) error stop 'rank3 mixed lane'
        end do
    end do
    ps15=-10.0_dp
    xr15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)=x(:,1)
    xr15(1,1,1,1,1,1,1,1,1,1,1,1,1,2,:)=x(:,2)
    call periodogram_into(xr15,ps15,f,axis=-1,status=st)
    if(st/=spectral_status_success) error stop 'rank15 periodogram status'
    if(maxval(abs(ps15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:) - &
        periodogram_power(x(:,1))))>tol) error stop 'rank15 periodogram mismatch'
    c15=cmplx(xr15,0.2_dp*xr15,dp)
    cp15=(-50.0_dp,3.0_dp)
    call csd_into(c15,c15,cp15,f,axis=-1,nperseg=4,noverlap=2,nfft=8,status=st)
    if(st/=spectral_status_success) error stop 'rank15 csd status'
    c=csd(cmplx(x(:,1),0.2_dp*x(:,1),dp),cmplx(x(:,1),0.2_dp*x(:,1),dp), &
        nperseg=4,noverlap=2,nfft=8)
    if(maxval(abs(cp15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-c%power))>tol) &
        error stop 'rank15 csd mismatch'
    padded_y=y(1:4,:)
    cross=(-77.0_dp,-33.0_dp)
    call csd_into(x,padded_y,cross(1:5,:),f,axis=1,nperseg=4,noverlap=2,nfft=8,status=st)
    if(st/=spectral_status_success) error stop 'unequal length status'
    do j=1,3
        c=csd(x(:,j),padded_y(:,j),nperseg=4,noverlap=2,nfft=8)
        if(maxval(abs(c%power-cross(1:5,j)))>tol) error stop 'unequal length value'
    end do
    strided(:,1:6:2)=x
    strided(:,2:6:2)=y
    strided_output=-11.0_dp
    call periodogram_into(strided(:,1:6:2),strided_output(:,2:6:2),f,axis=1,status=st)
    if(st/=spectral_status_success) error stop 'strided section status'
    do j=1,3
        p=periodogram(x(:,j))
        if(maxval(abs(p%power-strided_output(:,2*j)))>tol) error stop 'strided section value'
    end do
    if(any(strided_output(:,1:6:2)/=-11.0_dp)) error stop 'strided section untouched'
    short=10.0_dp
    call periodogram_into(x,short,f,axis=1,status=st)
    if(st/=spectral_status_invalid_input .or. any(short/=10.0_dp)) error stop 'shape unchanged'
    wrong=11.0_dp
    call periodogram_into(x,wrong,f,axis=0,status=st)
    if(st/=spectral_status_invalid_input .or. any(wrong/=11.0_dp)) error stop 'axis unchanged'
    expected=14.0_dp
    call welch_into(x,expected,f,axis=1,nperseg=4,noverlap=4,status=st)
    if(st/=spectral_status_invalid_input .or. any(expected/=14.0_dp)) error stop 'failure unchanged'
    print '(a)', 'test_signal_spectral_nd passed'
contains
    function periodogram_power(v) result(a)
        real(dp), intent(in) :: v(:) !! input signal
        real(dp), allocatable :: a(:)
        type(spectral_result) :: q
        q=periodogram(v)
        a=q%power
    end function periodogram_power
end program test_signal_spectral_nd
