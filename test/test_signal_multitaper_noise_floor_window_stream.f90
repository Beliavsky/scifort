! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_noise_floor_window_stream
    use scifort_kinds, only : dp
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal_spectral, only : spectral_status_success, spectral_status_invalid_input
    use scifort_signal_multitaper, only : multitaper_noise_floor_spectrogram, &
        multitaper_noise_floor_spectrogram_result, multitaper_noise_floor_spectrogram_channels_result, &
        multitaper_spectrogram, multitaper_spectrogram_channels_result, multitaper_noise_floor_channels_result, &
        multitaper_stream, multitaper_matrix_result
    implicit none
    integer, parameter :: n = 137, nc = 3, seg = 32, overlap = 13, nf = 45, width = 7
    integer, parameter :: chunks(9) = [0, 1, 19, 9, 31, 0, 7, 46, 24]
    real(dp) :: x(n,nc), t, pi, bad(2,nc)
    complex(dp) :: z(n,nc)
    type(multitaper_noise_floor_spectrogram_channels_result) :: batch, part, trans, empty
    type(multitaper_noise_floor_channels_result) :: mean_floor
    type(multitaper_noise_floor_spectrogram_result) :: single
    type(multitaper_spectrogram_channels_result) :: psd
    type(multitaper_stream) :: stream
    type(multitaper_matrix_result) :: avg
    real(dp), allocatable :: power(:,:,:), background(:,:,:), whitened(:,:,:), times(:)
    integer :: i, c, m, w, code, first, last, nframe, used
    real(dp), parameter :: tol = 5.0e-10_dp

    pi = acos(-1.0_dp)
    do c = 1, nc
        do i = 1, n
            t = real(i-1,dp)
            x(i,c) = 0.21_dp*t + 0.7_dp*sin(2*pi*(c+3)*t/31.0_dp) &
                + 0.6_dp*cos(2*pi*t/11.0_dp) + 0.3_dp*sin(0.17_dp*t*c)
            z(i,c) = cmplx(x(i,c),0.24_dp*t+0.4_dp*cos(0.33_dp*t*c),dp)
        end do
    end do
    do m = 1, 2
        do w = 1, 3
            select case (w)
            case (1)
                if (m == 1) then
                    batch = multitaper_noise_floor_spectrogram(x,3.0_dp,width,seg, &
                        noverlap=overlap,kmax=5,nfft=nf,fs=14.0_dp,weighting='equal',detrend_mode='linear')
                else
                    batch = multitaper_noise_floor_spectrogram(z,3.0_dp,width,seg, &
                        noverlap=overlap,kmax=5,nfft=nf,fs=14.0_dp,weighting='equal',detrend_mode='linear')
                end if
            case (2)
                if (m == 1) then
                    batch = multitaper_noise_floor_spectrogram(x,3.0_dp,width,seg, &
                        noverlap=overlap,kmax=5,nfft=nf,fs=14.0_dp,weighting='eigen',detrend_mode='linear')
                else
                    batch = multitaper_noise_floor_spectrogram(z,3.0_dp,width,seg, &
                        noverlap=overlap,kmax=5,nfft=nf,fs=14.0_dp,weighting='eigen',detrend_mode='linear')
                end if
            case default
                if (m == 1) then
                    batch = multitaper_noise_floor_spectrogram(x,3.0_dp,width,seg, &
                        noverlap=overlap,kmax=5,nfft=nf,fs=14.0_dp,weighting='adaptive',detrend_mode='linear')
                else
                    batch = multitaper_noise_floor_spectrogram(z,3.0_dp,width,seg, &
                        noverlap=overlap,kmax=5,nfft=nf,fs=14.0_dp,weighting='adaptive',detrend_mode='linear')
                end if
            end select
            if (batch%status /= spectral_status_success) error stop 'batch status'
            if (size(batch%power,2) /= 1+(n-seg)/(seg-overlap)) error stop 'batch frame count'
            if (size(batch%power,3) /= nc) error stop 'batch channel count'
            if (any(batch%background < 0.0_dp) .or. any(batch%whitened < 0.0_dp)) &
                error stop 'negative background/ratio'
            if (m == 1) then
                psd = multitaper_spectrogram(x,3.0_dp,seg,noverlap=overlap,kmax=5,nfft=nf, &
                    fs=14.0_dp,weighting=weighting_for(w),detrend_mode='linear')
                trans = multitaper_noise_floor_spectrogram(transpose(x),3.0_dp,width,seg, &
                    noverlap=overlap,kmax=5,nfft=nf,fs=14.0_dp, &
                    weighting=weighting_for(w),axis=2,detrend_mode='linear')
            else
                psd = multitaper_spectrogram(z,3.0_dp,seg,noverlap=overlap,kmax=5,nfft=nf, &
                    fs=14.0_dp,weighting=weighting_for(w),detrend_mode='linear')
                trans = multitaper_noise_floor_spectrogram(transpose(z),3.0_dp,width,seg, &
                    noverlap=overlap,kmax=5,nfft=nf,fs=14.0_dp, &
                    weighting=weighting_for(w),axis=2,detrend_mode='linear')
            end if
            if (psd%status /= spectral_status_success) error stop 'batch PSD status'
            if (maxval(abs(batch%power-psd%power)) > tol) error stop 'batch PSD disagreement'
            if (maxval(abs(batch%background-trans%background)) > tol) error stop 'transposed background'
            if (maxval(abs(batch%whitened-trans%whitened)) > tol) error stop 'transposed ratio'

            call stream%initialize(3.0_dp,seg,nc,code,noverlap=overlap,kmax=5, &
                fs=14.0_dp,nfft=nf,weighting=weighting_for(w),complex_input=m==2, &
                detrend_mode='linear')
            if (code /= spectral_status_success) error stop 'stream init'
            nframe = size(batch%power,2)
            allocate(power(size(batch%power,1),nframe,nc),background(size(batch%power,1),nframe,nc), &
                whitened(size(batch%power,1),nframe,nc),times(nframe))
            used = 0
            first = 1
            do c = 1, size(chunks)
                last = first+chunks(c)-1
                if (m == 1) then
                    if (mod(c,2) == 0) then
                        part = stream%push_noise_floor(transpose(x(first:last,:)),width,axis=2)
                    else
                        part = stream%push_noise_floor(x(first:last,:),width)
                    end if
                else
                    if (mod(c,2) == 0) then
                        part = stream%push_noise_floor(transpose(z(first:last,:)),width,axis=2)
                    else
                        part = stream%push_noise_floor(z(first:last,:),width)
                    end if
                end if
                if (part%status /= spectral_status_success) error stop 'stream push status'
                i = size(part%power,2)
                if (i > 0) then
                    power(:,used+1:used+i,:) = part%power
                    background(:,used+1:used+i,:) = part%background
                    whitened(:,used+1:used+i,:) = part%whitened
                    times(used+1:used+i) = part%times
                end if
                used = used+i
                first = last+1
            end do
            if (first /= n+1 .or. used /= nframe) error stop 'stream frame count'
            if (maxval(abs(power-batch%power)) > tol) error stop 'stream power'
            if (maxval(abs(background-batch%background)) > tol) error stop 'stream background'
            if (maxval(abs(whitened-batch%whitened)) > tol) error stop 'stream ratio'
            if (maxval(abs(times-batch%times)) > tol) error stop 'stream times'
            if (stream%sample_count() /= n) error stop 'stream count'
            avg = stream%mean_matrix()
            mean_floor = stream%mean_noise_floor(width)
            if (mean_floor%status /= spectral_status_success) error stop 'mean status'
            do i = 1, nc
                if (maxval(abs(mean_floor%power(:,i)-real(avg%power(:,i,i),dp))) > tol) &
                    error stop 'mean PSD'
                if (maxval(abs(mean_floor%power(:,i)-sum(batch%power(:,:,i),dim=2)/real(nframe,dp))) > tol) &
                    error stop 'mean batch PSD'
            end do
            empty = stream%push_noise_floor(x(1:0,:),4)
            if (empty%status /= spectral_status_invalid_input .or. &
                size(empty%power) /= 0) error stop 'even width accepted'
            if (stream%sample_count() /= n) error stop 'invalid width mutated state'
            bad = 0.0_dp
            bad(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
            if (m == 1) then
                empty = stream%push_noise_floor(bad,width)
            else
                empty = stream%push_noise_floor(cmplx(bad,0.0_dp,dp),width)
            end if
            if (empty%status /= spectral_status_invalid_input) error stop 'nonfinite accepted'
            if (stream%sample_count() /= n) error stop 'nonfinite mutated state'
            call stream%reset()
            if (stream%sample_count() /= 0 .or. stream%window_count() /= 0) error stop 'reset'
            mean_floor = stream%mean_noise_floor(width)
            if (mean_floor%status /= spectral_status_invalid_input) error stop 'empty mean accepted'
            deallocate(power,background,whitened,times)
        end do
    end do
    single = multitaper_noise_floor_spectrogram(x(:,1),3.0_dp,5,seg,nfft=nf)
    if (single%status /= spectral_status_success) error stop 'real single output'
    if (size(single%power,1) /= nf/2+1) error stop 'real single bins'
    single = multitaper_noise_floor_spectrogram(z(:,1),3.0_dp,5,seg,nfft=nf)
    if (single%status /= spectral_status_success) error stop 'complex single output'
    if (size(single%power,1) /= nf) error stop 'complex single bins'
    single = multitaper_noise_floor_spectrogram(x(:,1),3.0_dp,4,seg,nfft=nf)
    if (single%status /= spectral_status_invalid_input .or. size(single%power) /= 0) &
        error stop 'even single width accepted'
    single = multitaper_noise_floor_spectrogram(x(:,1),3.0_dp,99,seg,nfft=nf)
    if (single%status /= spectral_status_invalid_input) error stop 'wide single width accepted'
    ! Single-channel streaming uses the channels result even for rank-one input.
    ! Include even FFT length, noncontiguous samples, and a rejected axis.
    batch = multitaper_noise_floor_spectrogram(z(:,1:1),2.5_dp,5,16, &
        noverlap=8,kmax=4,fs=8.0_dp,nfft=24,weighting='eigen')
    call stream%initialize(2.5_dp,16,1,code,noverlap=8,kmax=4,fs=8.0_dp, &
        nfft=24,weighting='eigen',complex_input=.true.)
    if (code /= spectral_status_success) error stop 'single stream init'
    part = stream%push_noise_floor(z(:,1),5)
    if (part%status /= spectral_status_success) error stop 'single stream push'
    if (maxval(abs(part%power-batch%power)) > tol) error stop 'single stream power'
    if (maxval(abs(part%whitened-batch%whitened)) > tol) error stop 'single stream ratio'
    empty = stream%push_noise_floor(transpose(z(1:10,1:1)),5,axis=3)
    if (empty%status /= spectral_status_invalid_input) error stop 'invalid channel axis'
    if (stream%sample_count() /= n) error stop 'invalid axis changed stream count'
    call stream%reset()
    batch = multitaper_noise_floor_spectrogram(z(1::2,1:1),2.5_dp,5,16, &
        noverlap=8,kmax=4,fs=8.0_dp,nfft=24,weighting='eigen')
    part = stream%push_noise_floor(z(1::2,1),5)
    if (part%status /= spectral_status_success) error stop 'strided stream push'
    if (maxval(abs(part%background-batch%background)) > tol) error stop 'strided background'

    trans = multitaper_noise_floor_spectrogram(x,3.0_dp,5,seg,axis=3)
    if (trans%status /= spectral_status_invalid_input) error stop 'invalid batch axis'
    trans = multitaper_noise_floor_spectrogram(0.0_dp*x,3.0_dp,5,seg)
    if (trans%status /= spectral_status_success .or. any(trans%whitened /= 0.0_dp)) &
        error stop 'zero spectrum convention'
    print '(a)', 'test_signal_multitaper_noise_floor_window_stream passed'
contains
    function weighting_for(index) result(name)
        integer, intent(in) :: index
        character(8) :: name
        select case(index)
        case (1)
            name = 'equal'
        case (2)
            name = 'eigen'
        case default
            name = 'adaptive'
        end select
    end function weighting_for
end program test_signal_multitaper_noise_floor_window_stream
