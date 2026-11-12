! SPDX-License-Identifier: MIT
! Parameterized legacy ISTFT and COLA/NOLA checks using the public modules.
program test_signal_istft_parameter_windows
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : stft, get_window, stft_result, &
        spectral_status_success
    use scifort_signal_spectral_complex, only : stft_complex
    use scifort_signal_istft, only : istft, check_cola, check_nola, &
        istft_status_success, istft_status_invalid_input
    use scifort_signal_timefreq_nd, only : stft_into
    use scifort_signal_istft_nd, only : istft_into
    implicit none
    real(dp), parameter :: tol = 4.0e-11_dp
    real(dp) :: x(24), xr(24,2), output(24,2), wide_output(24,4)
    real(dp) :: xr14(1,1,1,1,1,1,1,1,1,1,1,1,1,24)
    real(dp) :: out14(1,1,1,1,1,1,1,1,1,1,1,1,1,24)
    real(dp) :: coeff(3), p(2), q(2)
    real(dp), allocatable :: y(:), times(:), freq(:), w(:)
    complex(dp) :: xc(24), zreal(6,2,13), zcomplex(10,2,13), &
        zw(6,4,13)
    complex(dp) :: yc(24,2)
    complex(dp) :: z15(1,1,1,1,1,1,1,1,1,1,1,1,1,6,13)
    complex(dp), allocatable :: reconstructed(:)
    type(stft_result) :: s
    integer :: i, j, k, status
    character(len=24) :: window

    x = [(sin(0.23_dp*real(i, dp)) + 0.14_dp*cos(0.91_dp*real(i, dp)) + &
        0.01_dp*real(i, dp), i=1,24)]
    xr(:,1) = x
    xr(:,2) = x*0.7_dp + 0.15_dp
    xc = cmplx(x, 0.4_dp*x - 0.12_dp, dp)
    coeff = [0.48_dp, 0.36_dp, 0.16_dp]
    p = [2.6_dp, 1.7_dp]
    q = [5.0_dp, 43.0_dp]

    do k=1,7
        select case(k)
        case(1)
            window='gaussian'
            s = stft(x, nperseg=8, noverlap=6, nfft=10, window=trim(window), &
                window_parameter=2.6_dp, boundary='even', scaling='psd')
            call istft(s%values, y, status, nperseg=8, noverlap=6, nfft=10, &
                window=trim(window), window_parameter=2.6_dp, scaling='psd')
        case(2)
            window='kaiser'
            s = stft(x, nperseg=8, noverlap=6, nfft=10, window=trim(window), &
                window_parameter=5.4_dp, boundary='even')
            call istft(s%values, y, status, nperseg=8, noverlap=6, nfft=10, &
                window=trim(window), window_parameter=5.4_dp)
        case(3)
            window='general_gaussian'
            s = stft(x, nperseg=8, noverlap=6, nfft=10, window=trim(window), &
                window_parameter=p(1), window_parameter2=p(2), boundary='even')
            call istft(s%values, y, status, nperseg=8, noverlap=6, nfft=10, &
                window=trim(window), window_parameter=p(1), window_parameter2=p(2))
        case(4)
            window='taylor'
            s = stft(x, nperseg=8, noverlap=6, nfft=10, window=trim(window), &
                window_parameter=q(1), window_parameter2=q(2), boundary='even')
            call istft(s%values, y, status, nperseg=8, noverlap=6, nfft=10, &
                window=trim(window), window_parameter=q(1), window_parameter2=q(2))
        case(5)
            window='general_cosine'
            s = stft(x, nperseg=8, noverlap=6, nfft=10, window=trim(window), &
                window_coefficients=coeff, boundary='even')
            call istft(s%values, y, status, nperseg=8, noverlap=6, nfft=10, &
                window=trim(window), window_coefficients=coeff)
        case(6)
            window='dpss'
            s = stft(x, nperseg=8, noverlap=6, nfft=10, window=trim(window), &
                window_parameter=2.1_dp, boundary='even')
            call istft(s%values, y, status, nperseg=8, noverlap=6, nfft=10, &
                window=trim(window), window_parameter=2.1_dp)
        case(7)
            window='hann'
            s = stft(x, nperseg=8, noverlap=6, nfft=10, window=trim(window), &
                boundary='even')
            call istft(s%values, y, status, nperseg=8, noverlap=6, nfft=10, window=trim(window))
        end select
        if (s%status /= spectral_status_success) error stop 'forward STFT window'
        if (status /= istft_status_success) error stop 'parameterized ISTFT'
        if (size(y) /= size(x)) error stop 'inverse output shape'
        if (maxval(abs(y-x)) > tol) error stop 'real signal round trip'
    end do

    s = stft_complex(xc, nperseg=8, noverlap=6, nfft=10, &
        window='gaussian', window_parameter=2.4_dp, boundary='odd', scaling='psd')
    if (s%status /= spectral_status_success) error stop 'complex forward STFT'
    call istft(s%values, reconstructed, status, nperseg=8, noverlap=6, &
        nfft=10, window='gaussian', window_parameter=2.4_dp, scaling='psd')
    if (status /= istft_status_success) error stop 'complex inverse STFT'
    if (size(reconstructed) /= 24) error stop 'complex output length'
    if (maxval(abs(reconstructed-xc)) > tol) error stop 'complex round trip'

    w = get_window('kaiser',8,parameter=5.4_dp)
    if (check_cola('kaiser',8,6,window_parameter=5.4_dp) .neqv. &
        check_cola(w,8,6)) error stop 'COLA named versus explicit'
    if (check_nola('kaiser',8,6,window_parameter=5.4_dp) .neqv. &
        check_nola(w,8,6)) error stop 'NOLA named versus explicit'
    if (check_nola('hann',8,4,window_parameter=4.0_dp,status=status)) &
        error stop 'invalid named window COLA accepted'
    if (status /= istft_status_invalid_input) error stop 'invalid NOLA status'

    call istft(s%values, reconstructed, status, window='hann', nperseg=8, &
        noverlap=6,nfft=10,window_parameter=3.0_dp)
    if (status /= istft_status_invalid_input) error stop 'unknown Hann parameter'
    call istft(s%values, reconstructed, status, window='gaussian',nperseg=8, &
        noverlap=6,nfft=10,window_parameter=3.0_dp,window_parameter2=2.0_dp)
    if (status /= istft_status_invalid_input) error stop 'unused second window parameter'
    call istft(s%values, reconstructed, status, window='gaussian',nperseg=8, &
        noverlap=6,nfft=10,window_parameter=3.0_dp,window_values=w)
    if (status /= istft_status_invalid_input) error stop 'conflicting window forms'
    call istft(s%values, reconstructed, status, window='general_cosine',nperseg=8, &
        noverlap=6,nfft=10)
    if (status /= istft_status_invalid_input) error stop 'missing cosine coefficients'

    call stft_into(xr,zreal,freq,times,window='gaussian',window_parameter=2.6_dp, &
        nperseg=8,noverlap=6,nfft=10,axis=1,boundary='constant',status=status)
    if (status /= spectral_status_success) error stop 'multi forward real'
    output=-333.0_dp
    call istft_into(zreal,output,times,window='gaussian',window_parameter=2.6_dp, &
        nperseg=8,noverlap=6,nfft=10,axis=1,status=status)
    if (status /= istft_status_success) error stop 'multi inverse real'
    if (maxval(abs(output-xr)) > tol) error stop 'multi real round trip'

    do j=1,2
        xc = cmplx(xr(:,j),0.34_dp*x,dp)
        s = stft_complex(xc,window='general_gaussian',window_parameter=1.8_dp, &
            window_parameter2=2.8_dp,nperseg=8,noverlap=6,nfft=10,boundary='odd')
        if (s%status/=spectral_status_success) error stop 'multi complex STFT'
        zcomplex(:,j,:) = s%values
    end do
    yc=(400.0_dp,-400.0_dp)
    call istft_into(zcomplex,yc,times,window='general_gaussian', &
        window_parameter=1.8_dp,window_parameter2=2.8_dp,nperseg=8, &
        noverlap=6,nfft=10,axis=1,status=status)
    if (status /= istft_status_success) error stop 'multi inverse complex'
    if (maxval(abs(real(yc,dp)-xr)) > tol) error stop 'multi complex real part'
    do j=1,2
        if (maxval(abs(aimag(yc(:,j))-0.34_dp*x))>tol) &
            error stop 'multi complex imaginary part'
    end do

    ! Failure must not modify a strided output, nor unused columns.
    zw=(0.0_dp,0.0_dp)
    zw(:,1:4:2,:)=zreal
    wide_output=9.0_dp
    call istft_into(zw(:,1:4:2,:),wide_output(:,1:4:2),times, &
        window='kaiser',window_parameter=2.0_dp,window_parameter2=3.0_dp, &
        nperseg=8,noverlap=6,nfft=10,axis=1,status=status)
    if (status /= istft_status_invalid_input) error stop 'multi invalid parameter'
    if (any(wide_output /= 9.0_dp)) error stop 'invalid modified output'
    if (size(times) /= 0) error stop 'invalid returned times'
    call istft_into(zw(:,1:4:2,:),wide_output(:,1:4:2),times, &
        window='gaussian',window_parameter=2.6_dp,nperseg=8,noverlap=6,nfft=10, &
        axis=1,status=status)
    if (status /= istft_status_success) error stop 'strided inverse'
    if (maxval(abs(wide_output(:,1:4:2)-xr))>tol) error stop 'strided inverse values'
    if (any(wide_output(:,2:4:2)/=9.0_dp)) error stop 'strided changed sentinel'

    xr14(1,1,1,1,1,1,1,1,1,1,1,1,1,:) = x
    call stft_into(xr14,z15,freq,times,window='gaussian',window_parameter=2.6_dp, &
        nperseg=8,noverlap=6,nfft=10,axis=-1,boundary='even',status=status)
    if (status /= spectral_status_success) error stop 'rank14 forward'
    out14=-999.0_dp
    call istft_into(z15,out14,times,window='gaussian',window_parameter=2.6_dp, &
        nperseg=8,noverlap=6,nfft=10,axis=-1,status=status)
    if (status /= istft_status_success) error stop 'rank14 inverse'
    if (maxval(abs(out14(1,1,1,1,1,1,1,1,1,1,1,1,1,:)-x))>tol) &
        error stop 'rank14 round trip'

    print '(a)', 'test_signal_istft_parameter_windows passed'
end program test_signal_istft_parameter_windows
