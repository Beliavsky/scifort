! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_istft_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal, only : dp, stft_into, istft_into, stft, &
        spectral_status_success, istft_status_success, istft_status_invalid_input, &
        istft_status_noninvertible, stft_result
    implicit none
    real(dp), parameter :: tol=5.0e-12_dp
    real(dp) :: x(16,3), xr(3,16), r(16,3), yr(3,16), yrbad(15,3)
    real(dp) :: x14(1,1,1,1,1,1,1,1,1,1,1,1,2,16)
    real(dp) :: y14(1,1,1,1,1,1,1,1,1,1,1,1,2,16)
    real(dp) :: x14first(1,1,1,1,1,1,1,1,1,1,1,1,16,2)
    real(dp) :: y14first(1,1,1,1,1,1,1,1,1,1,1,1,16,2)
    real(dp) :: strided(16,6)
    complex(dp) :: xc(3,16), yc(3,16), complex_from_real(16,3)
    complex(dp) :: z(5,3,5), zs(3,5,5), zc(3,8,5), zbad(5,3,5)
    complex(dp) :: z15(1,1,1,1,1,1,1,1,1,1,1,1,2,5,5)
    complex(dp) :: z15first(1,1,1,1,1,1,1,1,1,1,1,1,5,2,5)
    complex(dp) :: zstrided(5,6,5)
    complex(dp) :: z6(6,3,5), z_tf(5,6,3), z_ctf(3,5,6)
    complex(dp) :: z_tcf(5,3,6), z_ftc(6,5,3), z_cft(3,6,5)
    complex(dp) :: zc_tcf(5,3,8), z_padded(6,5,6)
    complex(dp) :: z15_shift(5,1,1,1,1,1,1,1,1,1,1,1,1,2,5)
    complex(dp), allocatable :: z1(:,:)
    real(dp), allocatable :: freq(:),t(:), inverse_times(:)
    type(stft_result) :: s
    integer :: i,j,b,f,st

    do j=1,3
        do i=1,16
            x(i,j)=sin(0.13_dp*i+0.21_dp*j)+0.22_dp*cos(0.37_dp*i-0.6_dp*j)
        end do
    end do
    r=-900.0_dp
    call stft_into(x,z,freq,t,fs=11.0_dp,nperseg=8,noverlap=4,nfft=8, &
        axis=1,status=st)
    if (st /= spectral_status_success) error stop 'forward rank2'
    call istft_into(z,r,inverse_times,fs=11.0_dp,nperseg=8,noverlap=4, &
        nfft=8,axis=1,status=st)
    if (st /= istft_status_success) error stop 'inverse rank3'
    if (maxval(abs(r-x))>tol) error stop 'real rank2 roundtrip'
    if (size(inverse_times)/=16) error stop 'time shape'
    if (abs(inverse_times(16)-15.0_dp/11.0_dp)>tol) error stop 'time grid'

    do j=1,3
        xr(j,:)=x(:,j)
        xc(j,:)=cmplx(x(:,j),0.33_dp*x(:,mod(j,3)+1),dp)
    end do
    call stft_into(xr,zs,freq,t,nperseg=8,noverlap=4,nfft=8,axis=-1,status=st)
    if (st /= spectral_status_success) error stop 'stft other axis'
    yr=-40.0_dp
    call istft_into(zs,yr,inverse_times,nperseg=8,noverlap=4,nfft=8, &
        axis=-1,status=st)
    if (st /= istft_status_success) error stop 'inverse other axis'
    if (maxval(abs(yr-xr))>tol) error stop 'inverse negative axis'
    call stft_into(xc,zc,freq,t,nperseg=8,noverlap=4,nfft=8,axis=2,status=st)
    if (st /= spectral_status_success) error stop 'complex forward'
    yc=(88.0_dp,24.0_dp)
    call istft_into(zc,yc,inverse_times,nperseg=8,noverlap=4,nfft=8, &
        axis=2,status=st)
    if (st /= istft_status_success) error stop 'complex inverse'
    if (maxval(abs(yc-xc))>tol) error stop 'complex inverse values'

    complex_from_real=(-80.0_dp,14.0_dp)
    call istft_into(z,complex_from_real,inverse_times,nperseg=8,noverlap=4, &
        nfft=8,axis=1,input_onesided=.true.,status=st)
    if (st /= istft_status_success) error stop 'onesided complex output'
    if (maxval(abs(complex_from_real-cmplx(x,0.0_dp,dp)))>tol) &
        error stop 'onesided complex values'

    ! Test maximum input rank with the frequency dimension next to time.
    x14=0.0_dp
    do j=1,2
        x14(1,1,1,1,1,1,1,1,1,1,1,1,j,:)=x(:,j)
    end do
    call stft_into(x14,z15,freq,t,nperseg=8,noverlap=4,nfft=8, &
        axis=-1,status=st)
    if (st /= spectral_status_success) error stop 'rank14 forward'
    y14=-100.0_dp
    call istft_into(z15,y14,inverse_times,nperseg=8,noverlap=4, &
        nfft=8,axis=-1,status=st)
    if (st /= istft_status_success) error stop 'rank15 inverse status'
    if (maxval(abs(y14-x14))>tol) error stop 'rank15 inverse values'

    ! Maximum-rank STFT with a frequency dimension before the channel axis.
    do j=1,2
        x14first(1,1,1,1,1,1,1,1,1,1,1,1,:,j)=x(:,j)
    end do
    call stft_into(x14first,z15first,freq,t,nperseg=8,noverlap=4, &
        nfft=8,axis=13,status=st)
    if (st /= spectral_status_success) error stop 'rank14 earlier axis forward'
    y14first=-42.0_dp
    call istft_into(z15first,y14first,inverse_times,nperseg=8,noverlap=4, &
        nfft=8,axis=13,status=st)
    if (st /= istft_status_success) error stop 'rank15 earlier axis inverse'
    if (maxval(abs(y14first-x14first))>tol) error stop 'rank15 earlier axis mismatch'

    ! Noncontiguous input/output, including unrelated destination elements.
    zstrided=(0.0_dp,0.0_dp)
    zstrided(:,2:6:2,:)=z
    strided=-999.0_dp
    call istft_into(zstrided(:,2:6:2,:),strided(:,1:6:2),inverse_times, &
        nperseg=8,noverlap=4,nfft=8,axis=1,status=st)
    if (st /= istft_status_success) error stop 'strided inverse'
    if (maxval(abs(strided(:,1:6:2)-x))>tol) error stop 'strided inverse values'
    if (any(strided(:,2:6:2)/=-999.0_dp)) error stop 'strided unrelated slots'

    r=-111.0_dp
    call istft_into(z,r,inverse_times,nperseg=8,noverlap=4,nfft=8,axis=0,status=st)
    if (st /= istft_status_invalid_input .or. any(r/=-111.0_dp)) &
        error stop 'invalid axis wrote result'
    if (size(inverse_times)/=0) error stop 'invalid axis times'

    ! A late failing channel must not commit earlier successful channels.
    zbad=z
    zbad(1,3,5)=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call istft_into(zbad,r,inverse_times,nperseg=8,noverlap=4,nfft=8, &
        axis=1,status=st)
    if (st /= istft_status_invalid_input .or. any(r/=-111.0_dp)) &
        error stop 'late channel failure wrote destination'
    if (size(inverse_times)/=0) error stop 'late channel failure returned times'
    call istft_into(z,r,inverse_times,nperseg=8,noverlap=4,nfft=7,axis=1,status=st)
    if (st /= istft_status_invalid_input .or. any(r/=-111.0_dp)) &
        error stop 'invalid fft wrote result'
    yrbad=-100.0_dp
    call istft_into(z,yrbad,inverse_times,nperseg=8,noverlap=4,nfft=8, &
        axis=1,status=st)
    if (st /= istft_status_invalid_input .or. any(yrbad/=-100.0_dp)) &
        error stop 'invalid output length wrote result'
    call istft_into(z,r,inverse_times,nperseg=8,noverlap=4,nfft=8, &
        axis=1,input_onesided=.false.,status=st)
    if (st /= istft_status_invalid_input .or. any(r/=-111.0_dp)) &
        error stop 'invalid real spectrum wrote result'

    ! Rank-two STFT input reconstructs a one-dimensional output.
    s=stft(x(:,1),window='hann',nperseg=8,noverlap=4,nfft=8)
    z1=s%values
    block
        real(dp) :: r1(16)
        r1=-19.0_dp
        call istft_into(z1,r1,inverse_times,window='hann',nperseg=8, &
            noverlap=4,nfft=8,status=st)
        if (st /= istft_status_success) error stop 'rank2 inverse'
        if (maxval(abs(r1-x(:,1)))>tol) error stop 'rank2 inverse values'
    end block

    block
        real(dp) :: r1(16)
        r1=-19.0_dp
        call istft_into(transpose(z1),r1,inverse_times,window='hann', &
            nperseg=8,noverlap=4,nfft=8,time_axis=1,freq_axis=2,status=st)
        if (st /= istft_status_success) error stop 'rank2 permuted axes'
        if (maxval(abs(r1-x(:,1)))>tol) error stop 'rank2 permuted axis values'
    end block

    s=stft(x(:,1),window='hann',nperseg=8,noverlap=4,nfft=8, &
        boundary='none' ,padded=.false.)
    if (s%status /= 0) error stop 'no boundary stft'
    z1=s%values
    block
        real(dp) :: r1(16)
        r1=13.0_dp
        call istft_into(z1,r1,inverse_times,window='hann',nperseg=8, &
            noverlap=4,nfft=8,boundary=.false.,status=st)
        if (st /= istft_status_noninvertible) error stop 'expected NOLA status'
        if (any(r1/=13.0_dp) .or. size(inverse_times)/=0) &
            error stop 'partial noninvertible reconstruction'
    end block

    ! Arbitrary pairs of input frequency/time axes (including reordered time
    ! first, middle, and last, and frequency after the time dimension).
    call stft_into(x,z6,freq,t,nperseg=8,noverlap=4,nfft=10,axis=1,status=st)
    if (st /= spectral_status_success) error stop 'frequency/time axis input'
    do j=1,3
        do b=1,5
            do f=1,6
                z_tf(b,f,j)=z6(f,j,b)
                z_ctf(j,b,f)=z6(f,j,b)
                z_tcf(b,j,f)=z6(f,j,b)
                z_ftc(f,b,j)=z6(f,j,b)
                z_cft(j,f,b)=z6(f,j,b)
            end do
        end do
    end do

    r=-70.0_dp
    call istft_into(z_tf,r,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=-3,freq_axis=-2,status=st)
    if (st/=istft_status_success .or. maxval(abs(r-x))>tol) &
        error stop 'time first frequency middle'
    r=-70.0_dp
    call istft_into(z_ftc,r,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=2,axis=1,status=st)
    if (st/=istft_status_success .or. maxval(abs(r-x))>tol) &
        error stop 'frequency first time middle'
    yr=-70.0_dp
    call istft_into(z_ctf,yr,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=2,freq_axis=3,status=st)
    if (st/=istft_status_success .or. maxval(abs(yr-xr))>tol) &
        error stop 'frequency last time middle'
    yr=-70.0_dp
    call istft_into(z_tcf,yr,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=1,freq_axis=3,axis=-1,status=st)
    if (st/=istft_status_success .or. maxval(abs(yr-xr))>tol) &
        error stop 'time first frequency last'
    yr=-70.0_dp
    call istft_into(z_cft,yr,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        freq_axis=2,status=st)
    if (st/=istft_status_success .or. maxval(abs(yr-xr))>tol) &
        error stop 'time last frequency middle'

    ! The same layout works with complex two-sided input and output.
    do j=1,3
        do b=1,5
            do f=1,8
                zc_tcf(b,j,f)=zc(j,f,b)
            end do
        end do
    end do
    yc=(-77.0_dp,11.0_dp)
    call istft_into(zc_tcf,yc,inverse_times,nperseg=8,noverlap=4, &
        nfft=8,time_axis=1,freq_axis=3,status=st)
    if (st/=istft_status_success .or. maxval(abs(yc-xc))>tol) &
        error stop 'complex time first frequency last'

    ! A rank-15 source can have its time axis first and frequency last.
    do j=1,2
        do f=1,5
            z15_shift(:,1,1,1,1,1,1,1,1,1,1,1,1,j,f) = &
                z15(1,1,1,1,1,1,1,1,1,1,1,1,j,f,:)
        end do
    end do
    y14=-70.0_dp
    call istft_into(z15_shift,y14,inverse_times,nperseg=8,noverlap=4, &
        nfft=8,time_axis=1,freq_axis=-1,status=st)
    if (st/=istft_status_success .or. maxval(abs(y14-x14))>tol) &
        error stop 'rank15 permuted transform axes'
    y14=-70.0_dp
    call istft_into(z15_shift,y14,inverse_times,nperseg=8,noverlap=4, &
        nfft=8,time_axis=1,status=st)
    if (st/=istft_status_success .or. maxval(abs(y14-x14))>tol) &
        error stop 'rank15 inferred last frequency axis' 

    ! Strided rank-three input and noncontiguous output.
    z_padded=(-19.0_dp,17.0_dp)
    z_padded(1:6:2,:,:)=z_ctf
    strided=-777.0_dp
    call istft_into(z_padded(1:6:2,:,:),strided(1:6:2,:),inverse_times, &
        nperseg=8,noverlap=4,nfft=10,time_axis=2,freq_axis=3,status=st)
    ! The output section has shape (3,6), not the required (3,16).
    if (st/=istft_status_invalid_input .or. any(strided/=-777.0_dp)) &
        error stop 'strided output shape not validated'
    strided=-777.0_dp
    call istft_into(z_padded(1:6:2,:,:),strided, inverse_times, &
        nperseg=8,noverlap=4,nfft=10,time_axis=2,freq_axis=3,status=st)
    if (st/=istft_status_invalid_input .or. any(strided/=-777.0_dp)) &
        error stop 'mismatched strided output wrote values'
    ! Use the matching channel-first strided signal slice.
    block
        real(dp) :: dest(6,16)
        dest=-777.0_dp
        call istft_into(z_padded(1:6:2,:,:),dest(1:6:2,:),inverse_times, &
            nperseg=8,noverlap=4,nfft=10,time_axis=2,freq_axis=3,status=st)
        if (st/=istft_status_success) error stop 'strided non-last axes'
        if (maxval(abs(dest(1:6:2,:)-xr))>tol) error stop 'strided non-last values'
        if (any(dest(2:6:2,:)/=-777.0_dp)) error stop 'strided other channels'
    end block

    r=-111.0_dp
    call istft_into(z_tf,r,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=1,freq_axis=1,status=st)
    if (st/=istft_status_invalid_input .or. any(r/=-111.0_dp) .or. size(inverse_times)/=0) &
        error stop 'duplicate transform axes not rejected'
    call istft_into(z_tf,r,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=0,freq_axis=2,status=st)
    if (st/=istft_status_invalid_input .or. any(r/=-111.0_dp)) &
        error stop 'bad time axis not rejected'
    call istft_into(z_tf,r,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=1,freq_axis=4,status=st)
    if (st/=istft_status_invalid_input .or. any(r/=-111.0_dp)) &
        error stop 'bad frequency axis not rejected'
    call istft_into(z_tf,r,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=1,freq_axis=2,axis=2,status=st)
    if (st/=istft_status_invalid_input .or. any(r/=-111.0_dp)) &
        error stop 'contradictory sample/frequency axes not rejected'
    call istft_into(z_ctf,r,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=2,freq_axis=3,status=st)
    if (st/=istft_status_invalid_input .or. any(r/=-111.0_dp)) &
        error stop 'wrong channel dimensions not rejected'
    ! A corrupt coefficient in the final channel must not commit an output
    ! reconstructed from earlier, otherwise valid, channels.
    z_tf(5,6,3)=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call istft_into(z_tf,r,inverse_times,nperseg=8,noverlap=4,nfft=10, &
        time_axis=1,freq_axis=2,status=st)
    if (st/=istft_status_invalid_input .or. any(r/=-111.0_dp) .or. size(inverse_times)/=0) &
        error stop 'late-channel failure with reordered axes wrote destination' 

    print '(a)', 'test_signal_istft_nd passed'
end program test_signal_istft_nd
