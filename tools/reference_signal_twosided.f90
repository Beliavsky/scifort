! SPDX-License-Identifier: MIT
! Fixed deterministic real and complex signals for independent SciPy comparisons.
program reference_signal_twosided
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : periodogram, welch, csd, stft, &
        spectrogram, spectral_result, cross_spectral_result, stft_result, &
        spectrogram_result, spectral_status_success
    use scifort_signal_spectral_complex, only : spectrogram_complex
    implicit none

    integer, parameter :: n = 37
    integer, parameter :: fft_lengths(4) = [8,9,12,17]
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: x(n), y(n)
    complex(dp) :: z(n)
    real(dp), allocatable :: fs(:), vs(:), times(:), vals(:,:), shape_param
    complex(dp), allocatable :: cs(:), cv(:,:)
    type(spectral_result) :: power
    type(cross_spectral_result) :: cross
    type(spectrogram_result) :: picture
    type(stft_result) :: transform
    integer :: k, fi, wi, si, oi, family, nf, count, frame
    logical :: onesided
    character(:), allocatable :: wname, scaling, trend

    do k = 1,n
        x(k) = 0.18_dp*real(k,dp)+0.7_dp*sin(0.48_dp*real(k,dp)) + &
            0.36_dp*cos(0.17_dp*real(k*k,dp))
        y(k) = -0.06_dp*real(k,dp)+0.8_dp*cos(0.26_dp*real(k,dp)) + &
            0.25_dp*sin(0.37_dp*real(k,dp))
        z(k) = cmplx(x(k), y(k), dp)
    end do
    count = 0
    do family = 1, 6
        do fi = 1, size(fft_lengths)
            nf = fft_lengths(fi)
            do wi = 1, 2
                wname = 'hann'
                if (allocated(shape_param)) deallocate(shape_param)
                if (wi == 2) then
                    wname = 'gaussian'
                    shape_param = 2.25_dp
                end if
                do si = 1, 2
                    scaling = 'density'
                    if (si == 2) scaling = 'spectrum'
                    do oi = 1, 2
                        onesided = oi == 1
                        trend = 'constant'
                        if (oi == 2) trend = 'linear'
                        select case (family)
                        case (1)
                            power = periodogram(x, fs=25.0_dp, window=wname, nfft=nf, &
                                detrend=trend, scaling=scaling, window_parameter=shape_param, &
                                return_onesided=onesided)
                            if (power%status /= spectral_status_success) error stop 'periodogram failed'
                            fs = power%frequencies
                            vs = power%power
                        case (2)
                            power = welch(x, fs=25.0_dp, window=wname, nperseg=8, noverlap=3, &
                                nfft=nf, detrend=trend, scaling=scaling, average=average_mode(wi), &
                                window_parameter=shape_param, return_onesided=onesided)
                            if (power%status /= spectral_status_success) error stop 'welch failed'
                            fs = power%frequencies
                            vs = power%power
                        case (3)
                            cross = csd(x, y, fs=25.0_dp, window=wname, nperseg=8, noverlap=3, &
                                nfft=nf, detrend=trend, scaling=scaling, average=average_mode(wi), &
                                window_parameter=shape_param, return_onesided=onesided)
                            if (cross%status /= spectral_status_success) error stop 'csd failed'
                            fs = cross%frequencies
                            cs = cross%power
                        case (4)
                            transform = stft(x, fs=25.0_dp, window=wname, nperseg=8, noverlap=3, &
                                nfft=nf, detrend=trend, boundary=stft_boundary(wi), padded=wi==2, &
                                scaling=stft_scale(si), window_parameter=shape_param, &
                                return_onesided=onesided)
                            if (transform%status /= spectral_status_success) error stop 'stft failed'
                            fs = transform%frequencies
                            times = transform%times
                            cv = transform%values
                        case (5)
                            picture = spectrogram(x, fs=25.0_dp, window=wname, nperseg=8, noverlap=3, &
                                nfft=nf, detrend=trend, scaling=scaling, window_parameter=shape_param, &
                                return_onesided=onesided)
                            if (picture%status /= spectral_status_success) error stop 'spectrogram failed'
                            fs = picture%frequencies
                            times = picture%times
                            vals = picture%power
                        case (6)
                            picture = spectrogram_complex(z, fs=25.0_dp, window=wname, nperseg=8, &
                                noverlap=3, nfft=nf, detrend=trend, scaling=scaling, &
                                window_parameter=shape_param, return_onesided=onesided)
                            if (picture%status /= spectral_status_success) error stop 'complex spectrogram failed'
                            fs = picture%frequencies
                            times = picture%times
                            vals = picture%power
                        end select
                        count = count+1
                        frame = 0
                        select case (family)
                        case (4)
                            frame = size(cv,2)
                        case (5,6)
                            frame = size(vals,2)
                        end select
                        print '(a,7(1x,i0))', 'CASE',family,nf,wi,si,oi,size(fs),frame
                        do k=1,size(fs)
                            print '(es26.17)',fs(k)
                        end do
                        select case (family)
                        case (1,2)
                            do k=1,size(vs)
                                print '(es26.17)',vs(k)
                            end do
                        case (3)
                            do k=1,size(cs)
                                print '(2(es26.17,1x))',real(cs(k),dp),aimag(cs(k))
                            end do
                        case (4)
                            do k=1,frame
                                print '(es26.17)',times(k)
                            end do
                            do frame=1,size(cv,2)
                                do k=1,size(cv,1)
                                    print '(2(es26.17,1x))',real(cv(k,frame),dp),aimag(cv(k,frame))
                                end do
                            end do
                        case (5,6)
                            do k=1,frame
                                print '(es26.17)',times(k)
                            end do
                            do frame=1,size(vals,2)
                                do k=1,size(vals,1)
                                    print '(es26.17)',vals(k,frame)
                                end do
                            end do
                        end select
                    end do
                end do
            end do
        end do
    end do
    if (count /= 192) error stop 'case count wrong'
contains
    function average_mode(w) result(mode)
        integer, intent(in) :: w !! median for second window
        character(:), allocatable :: mode
        mode='mean'
        if (w == 2) mode='median'
    end function average_mode

    function stft_boundary(w) result(mode)
        integer, intent(in) :: w !! boundary modes vary with window family
        character(:), allocatable :: mode
        mode='even'
        if (w == 2) mode='odd'
    end function stft_boundary

    function stft_scale(i) result(mode)
        integer, intent(in) :: i !! scaling selector
        character(:), allocatable :: mode
        mode='psd'
        if (i == 2) mode='spectrum'
    end function stft_scale
end program reference_signal_twosided
