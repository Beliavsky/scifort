! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_stft_boundaries
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : stft, stft_result, spectral_status_success, &
        spectral_status_invalid_input
    use scifort_signal_spectral_complex, only : stft_complex
    use scifort_signal_timefreq_nd, only : stft_into
    implicit none
    real(dp), parameter :: tol = 1.0e-12_dp
    real(dp) :: x(4), a(4,2), ar(4,4), fexp(3)
    complex(dp) :: z(4), z2(4,2), yz(4,2,3), ys(3,4,3)
    complex(dp) :: ybad(3,2,3), y14(1,1,1,1,1,1,1,1,1,1,1,1,1,3,3)
    real(dp) :: xr14(1,1,1,1,1,1,1,1,1,1,1,1,1,4)
    real(dp), allocatable :: f(:), times(:)
    type(stft_result) :: r, c
    integer :: st, i, j
    character(8), parameter :: modes(5) = [character(8) :: &
        'even', 'odd', 'constant', 'zeros', 'none']
    real(dp), parameter :: dc(5) = [2.0_dp, 0.5_dp, 1.25_dp, 0.75_dp, 2.5_dp]

    x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    z = cmplx(x, x*x, dp)
    fexp = [0.0_dp, 1.0_dp, 2.0_dp]
    do i = 1, 5
        r = stft(x, window='boxcar', nperseg=4, noverlap=2, nfft=4, &
            fs=2.0_dp, detrend='none', boundary=trim(modes(i)), padded=.false.)
        if (r%status /= spectral_status_success) error stop 'real boundary status'
        if (abs(real(r%values(1,1),dp)-dc(i)) > tol) error stop 'real boundary DC'
        if (abs(aimag(r%values(1,1))) > tol) error stop 'real DC imaginary'
        if (i < 5) then
            if (size(r%times) /= 3) error stop 'extended frame count'
            if (maxval(abs(r%times-fexp)) > tol) error stop 'extended time grid'
        else
            if (size(r%times) /= 1) error stop 'no-extension frame count'
            if (abs(r%times(1)-1.0_dp) > tol) error stop 'no-extension time'
        end if
        c = stft_complex(z, window='boxcar', nperseg=4, noverlap=2, nfft=4, &
            fs=2.0_dp, detrend='none', boundary=trim(modes(i)), padded=.false.)
        if (c%status /= spectral_status_success) error stop 'complex boundary status'
        if (maxval(abs(c%times-r%times)) > tol) error stop 'real/complex time mismatch'
        if (maxval(abs(c%values(1:3,:)-r%values)) < 0.01_dp) &
            error stop 'complex imaginary components lost'
        ! A complex signal with no imaginary part must reproduce the one-sided result.
        c = stft_complex(cmplx(x,0.0_dp,dp), window='boxcar', nperseg=4, &
            noverlap=2, nfft=4, fs=2.0_dp, detrend='none', &
            boundary=trim(modes(i)), padded=.false.)
        if (maxval(abs(c%values(1:3,:)-r%values)) > tol) &
            error stop 'real/complex shared extension'
    end do

    r = stft(x,boundary='wrap')
    if (r%status /= spectral_status_invalid_input) error stop 'unsupported boundary accepted'
    c = stft_complex(z,boundary='wrap')
    if (c%status /= spectral_status_invalid_input) error stop 'complex unsupported boundary'

    do j=1,2
        a(:,j) = x + real(j,dp)
        z2(:,j) = cmplx(a(:,j),-0.2_dp*a(:,j),dp)
    end do
    yz = (-999.0_dp,0.0_dp)
    call stft_into(z2,yz,f,times,window='boxcar',nperseg=4,noverlap=2,nfft=4, &
        boundary='odd',detrend='none',padded=.false.,axis=1,status=st)
    if (st /= spectral_status_success) error stop 'multichannel boundary status'
    do j=1,2
        c = stft_complex(z2(:,j),window='boxcar',nperseg=4,noverlap=2,nfft=4, &
            boundary='odd',detrend='none',padded=.false.)
        if (maxval(abs(yz(:,j,:)-c%values)) > tol) error stop 'multichannel boundary values'
    end do
    ! The same operation must support noncontiguous input and output sections.
    ar = -200.0_dp
    ar(:,1:4:2) = a
    ys = cmplx(-100.0_dp,5.0_dp,dp)
    call stft_into(ar(:,1:4:2),ys(:,1:4:2,:),f,times,window='boxcar', &
        nperseg=4,noverlap=2,nfft=4,boundary='even',detrend='none',padded=.false., &
        axis=-2,status=st)
    if (st /= spectral_status_success) error stop 'strided extension status'
    if (any(ys(:,2:4:2,:) /= cmplx(-100.0_dp,5.0_dp,dp))) &
        error stop 'strided extension unrelated channels'
    do j=1,2
        r = stft(a(:,j),window='boxcar',nperseg=4,noverlap=2,nfft=4, &
            boundary='even',detrend='none',padded=.false.)
        if (maxval(abs(ys(:,2*j-1,:)-r%values)) > tol) &
            error stop 'strided extension values'
    end do
    ybad=(9.0_dp,-5.0_dp)
    call stft_into(a,ybad,f,times,nperseg=4,noverlap=2,nfft=4, &
        boundary='bad',status=st)
    if (st /= spectral_status_invalid_input) error stop 'invalid wrapper boundary status'
    if (any(ybad /= cmplx(9.0_dp,-5.0_dp,dp))) &
        error stop 'invalid boundary modified caller output'
    if (size(f) /= 0 .or. size(times) /= 0) error stop 'invalid coordinate outputs'

    xr14(1,1,1,1,1,1,1,1,1,1,1,1,1,:)=x
    y14=(0.0_dp,0.0_dp)
    call stft_into(xr14,y14,f,times,window='boxcar',nperseg=4,noverlap=2,nfft=4, &
        detrend='none',boundary='constant',padded=.false.,axis=-1,status=st)
    if (st /= spectral_status_success) error stop 'rank14 boundary status'
    r=stft(x,window='boxcar',nperseg=4,noverlap=2,nfft=4, &
        detrend='none',boundary='constant',padded=.false.)
    if (maxval(abs(y14(1,1,1,1,1,1,1,1,1,1,1,1,1,:,:)-r%values)) > tol) &
        error stop 'rank14 boundary values'
    print '(a)', 'test_signal_stft_boundaries passed'
end program test_signal_stft_boundaries
