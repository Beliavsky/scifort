! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style complex, two-sided 1D spectral estimators. Independent code;
! standard Welch/cross-spectrum and Fourier-transform definitions.
module scifort_signal_spectral_complex
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, &
        stft_result, spectrogram_result, spectral_status_success, &
        spectral_status_invalid_input, spectral_setup, forward_fft, &
        median_bias, array_median, stft_edge_mode_valid, fill_stft_edges, spectral_unwrap_phase
    implicit none
    private

    public :: periodogram_complex, welch_complex, csd_complex, &
        coherence_complex, stft_complex, spectrogram_complex
    public :: csd_real_complex, csd_complex_real, &
        coherence_real_complex, coherence_complex_real
    public :: spectral_complex_detrend_callback

    abstract interface
        subroutine spectral_complex_detrend_callback(segment)
            import dp
            complex(dp), intent(inout) :: segment(:) !! one unwindowed complex segment
        end subroutine spectral_complex_detrend_callback
    end interface

contains

    pure logical function finite_complex(x) result(ok)
        complex(dp), intent(in) :: x(:) !! complex input signal
        ok = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function finite_complex

    subroutine detrend_complex_segment(seg, method, ok)
        complex(dp), intent(inout) :: seg(:) !! unwindowed segment
        character(*), intent(in) :: method !! 'constant', 'linear', or 'none'
        logical, intent(out) :: ok !! whether the selected mode is supported
        complex(dp) :: av, slope
        real(dp) :: middle, denom
        integer :: i, n

        n = size(seg)
        ok = .true.
        select case (trim(method))
        case ('none', 'false')
        case ('constant')
            if (n > 0) seg = seg - sum(seg)/real(n,dp)
        case ('linear')
            if (n == 0) return
            av = sum(seg)/real(n,dp)
            middle = real(n+1,dp)/2.0_dp
            slope = (0.0_dp,0.0_dp)
            if (n > 1) then
                denom = real(n,dp)*(real(n,dp)**2 - 1.0_dp)/12.0_dp
                do i = 1, n
                    slope = slope + (real(i,dp)-middle)*(seg(i)-av)
                end do
                slope = slope/denom
            end if
            do i = 1, n
                seg(i) = seg(i) - av - (real(i,dp)-middle)*slope
            end do
        case default
            ok = .false.
        end select
    end subroutine detrend_complex_segment

    subroutine transform_complex_segment(x, offset, nw, nf, w, trend, ft, ok, detrend_callback)
        complex(dp), intent(in) :: x(:) !! signal, including any padding
        integer, intent(in) :: offset !! one-based segment start
        integer, intent(in) :: nw !! window length
        integer, intent(in) :: nf !! FFT length
        real(dp), intent(in) :: w(:) !! real analysis window
        character(*), intent(in) :: trend !! detrend mode
        complex(dp), intent(out) :: ft(:) !! length-nf Fourier output
        logical, intent(out) :: ok !! valid segment and trend
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        complex(dp) :: local(nw)

        ok = .false.
        if (offset < 1 .or. offset+nw-1 > size(x)) return
        local = x(offset:offset+nw-1)
        if (present(detrend_callback)) then
            call detrend_callback(local)
            ok = finite_complex(local)
        else
            call detrend_complex_segment(local,trend,ok)
        end if
        if (.not. ok) return
        ft = (0.0_dp,0.0_dp)
        ft(1:nw) = local*w
        call forward_fft(ft)
    end subroutine transform_complex_segment

    subroutine two_sided_frequencies(n, fs, freqs)
        integer, intent(in) :: n !! transform length
        real(dp), intent(in) :: fs !! sample rate
        real(dp), intent(out) :: freqs(:) !! numpy.fft.fftfreq order
        integer :: k, signed_k

        do k = 0, n-1
            signed_k = k
            if (k > (n-1)/2) signed_k = k - n
            freqs(k+1) = real(signed_k,dp)*fs/real(n,dp)
        end do
    end subroutine two_sided_frequencies

    function csd_complex(x, y, fs, window, nperseg, noverlap, nfft, detrend, &
            scaling, average, return_onesided, &
            window_parameter, window_parameter2, window_coefficients, window_values, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! first complex signal
        complex(dp), intent(in) :: y(:) !! second complex signal
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! default 'hann'
        integer, intent(in), optional :: nperseg !! segment size
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: nfft !! FFT length
        character(*), intent(in), optional :: detrend !! 'constant', 'linear', 'none'
        character(*), intent(in), optional :: scaling !! 'density' or 'spectrum'
        character(*), intent(in), optional :: average !! 'mean' or 'median'
        logical, intent(in), optional :: return_onesided !! ignored for complex signals, as in SciPy

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(cross_spectral_result) :: r
        real(dp), allocatable :: w(:), original(:), imagpart(:)
        complex(dp), allocatable :: xx(:), yy(:), fx(:), fy(:), spectra(:,:)
        real(dp) :: rate, factor, correction
        integer :: n, nw, no, nf, segments, hop, j, k
        logical :: ok
        character(:), allocatable :: trend, avg

        allocate(r%frequencies(0), r%power(0))
        if (size(x) < 1 .or. size(y) < 1) then
            if (size(x) == 0 .and. size(y) == 0) r%status = spectral_status_success
            return
        end if
        if (.not. (finite_complex(x) .and. finite_complex(y))) return
        n = max(size(x),size(y))
        call spectral_setup(n, fs, window, nperseg, noverlap, nfft, scaling, &
            'hann', 0, nw, no, nf, rate, w, factor, ok, &
            window_parameter, window_parameter2, window_coefficients, window_values)
        if (.not. ok) return
        trend = 'constant'
        if (present(detrend)) trend = trim(detrend)
        if (present(detrend_callback) .and. present(detrend)) return
        if (.not. present(detrend_callback)) then
            if (trend /= 'constant' .and. trend /= 'linear' .and. &
                trend /= 'none' .and. trend /= 'false') return
        end if
        avg = 'mean'
        if (present(average)) avg = trim(average)
        if (avg /= 'mean' .and. avg /= 'median') return
        hop = nw - no
        segments = 1 + (n-nw)/hop
        allocate(xx(n),yy(n),fx(nf),fy(nf),spectra(nf,segments))
        xx = (0.0_dp,0.0_dp)
        yy = (0.0_dp,0.0_dp)
        xx(1:size(x)) = x
        yy(1:size(y)) = y
        do j = 1, segments
            call transform_complex_segment(xx,1+(j-1)*hop,nw,nf,w,trend,fx,ok,detrend_callback)
            if (.not. ok) return
            call transform_complex_segment(yy,1+(j-1)*hop,nw,nf,w,trend,fy,ok,detrend_callback)
            if (.not. ok) return
            spectra(:,j) = conjg(fx)*fy*factor
        end do
        deallocate(r%frequencies,r%power)
        allocate(r%frequencies(nf),r%power(nf))
        call two_sided_frequencies(nf,rate,r%frequencies)
        if (avg == 'mean') then
            r%power = sum(spectra,dim=2)/real(segments,dp)
        else
            allocate(original(segments),imagpart(segments))
            correction = median_bias(segments)
            do k = 1, nf
                original = real(spectra(k,:))
                imagpart = aimag(spectra(k,:))
                r%power(k) = cmplx(array_median(original),array_median(imagpart),dp)/correction
            end do
        end if
        r%status = spectral_status_success
    end function csd_complex

    function periodogram_complex(x, fs, window, nfft, detrend, scaling, return_onesided, &
            window_parameter, window_parameter2, window_coefficients, window_values, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! complex input
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! default boxcar
        integer, intent(in), optional :: nfft !! FFT size, can truncate input
        character(*), intent(in), optional :: detrend !! detrend mode
        character(*), intent(in), optional :: scaling !! density or spectrum
        logical, intent(in), optional :: return_onesided !! ignored for complex

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_result) :: r
        type(cross_spectral_result) :: cross
        character(:), allocatable :: win
        integer :: nuse

        allocate(r%frequencies(0),r%power(0))
        if (size(x) == 0) then
            r%status = spectral_status_success
            return
        end if
        nuse = size(x)
        if (present(nfft)) then
            if (nfft < 1) return
            nuse = min(nuse,nfft)
        end if
        if (.not. present(window_values)) win = 'boxcar'
        if (present(window)) win = trim(window)
        cross = csd_complex(x(:nuse),x(:nuse),fs=fs,window=win,nperseg=nuse, &
            noverlap=0,nfft=nfft,detrend=detrend,scaling=scaling, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
        r%frequencies = cross%frequencies
        r%power = real(cross%power,dp)
        r%status = cross%status
    end function periodogram_complex

    function welch_complex(x, fs, window, nperseg, noverlap, nfft, detrend, &
            scaling, average, return_onesided, &
            window_parameter, window_parameter2, window_coefficients, window_values, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! complex input
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! default hann
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap length
        integer, intent(in), optional :: nfft !! FFT length
        character(*), intent(in), optional :: detrend !! detrend mode
        character(*), intent(in), optional :: scaling !! density or spectrum
        character(*), intent(in), optional :: average !! mean or median
        logical, intent(in), optional :: return_onesided !! ignored for complex

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_result) :: r
        type(cross_spectral_result) :: cross

        cross = csd_complex(x,x,fs,window,nperseg,noverlap,nfft,detrend,scaling,average, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
        r%frequencies = cross%frequencies
        r%power = real(cross%power,dp)
        r%status = cross%status
    end function welch_complex

    function coherence_complex(x,y,fs,window,nperseg,noverlap,nfft,detrend, &
            window_parameter, window_parameter2, window_coefficients, window_values, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! first complex signal
        complex(dp), intent(in) :: y(:) !! second complex signal
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! window name
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap length
        integer, intent(in), optional :: nfft !! transform length
        character(*), intent(in), optional :: detrend !! detrending mode

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_result) :: r
        type(cross_spectral_result) :: c, pxx, pyy
        integer :: i
        real(dp) :: denom

        allocate(r%frequencies(0),r%power(0))
        c = csd_complex(x,y,fs,window,nperseg,noverlap,nfft,detrend, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
        if (c%status /= spectral_status_success) return
        pxx = csd_complex(x,x,fs,window,nperseg,noverlap,nfft,detrend, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
        pyy = csd_complex(y,y,fs,window,nperseg,noverlap,nfft,detrend, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
        if (pxx%status /= spectral_status_success .or. &
            pyy%status /= spectral_status_success) return
        r%frequencies = c%frequencies
        deallocate(r%power)
        allocate(r%power(size(c%power)))
        do i = 1, size(r%power)
            denom = real(pxx%power(i),dp)*real(pyy%power(i),dp)
            r%power(i) = ieee_value(0.0_dp,ieee_quiet_nan)
            if (denom > tiny(1.0_dp)) r%power(i) = abs(c%power(i))**2/denom
        end do
        r%status = spectral_status_success
    end function coherence_complex

    function stft_complex(x,fs,window,nperseg,noverlap,nfft,detrend,boundary, &
            padded,scaling,return_onesided, &
            window_parameter, window_parameter2, window_coefficients, window_values, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! complex input
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! default hann
        integer, intent(in), optional :: nperseg !! segment size
        integer, intent(in), optional :: noverlap !! overlap length
        integer, intent(in), optional :: nfft !! FFT size
        character(*), intent(in), optional :: detrend !! detrend mode
        character(*), intent(in), optional :: boundary !! zeros (default), even, odd, constant, none
        logical, intent(in), optional :: padded !! append zeros for last segment
        character(*), intent(in), optional :: scaling !! 'spectrum' or 'psd'
        logical, intent(in), optional :: return_onesided !! ignored for complex

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(stft_result) :: r

        r = compute_complex_stft(x,fs,window,nperseg,noverlap,nfft,detrend, &
            boundary,padded,scaling,.false., &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
    end function stft_complex

    function spectrogram_complex(x,fs,window,nperseg,noverlap,nfft,detrend,scaling, &
            window_parameter, window_parameter2, window_coefficients, window_values, &
            return_onesided, mode, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! complex signal
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! default Tukey
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap length
        integer, intent(in), optional :: nfft !! FFT length
        character(*), intent(in), optional :: detrend !! default constant
        character(*), intent(in), optional :: scaling !! density or spectrum

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        logical, intent(in), optional :: return_onesided !! complex signals always use two-sided output
        character(*), intent(in), optional :: mode !! psd (default), complex, magnitude, angle or phase
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectrogram_result) :: r
        type(stft_result) :: temp
        character(:), allocatable :: scale, chosen_mode

        chosen_mode = 'psd'
        if (present(mode)) chosen_mode = trim(mode)
        select case (chosen_mode)
        case ('psd', 'complex', 'magnitude', 'angle', 'phase')
        case default
            allocate(r%frequencies(0), r%times(0), r%power(0,0))
            return
        end select
        scale = 'density'
        if (present(scaling)) scale = trim(scaling)
        temp = compute_complex_stft(x,fs,window=window,nperseg=nperseg, &
            noverlap=noverlap,nfft=nfft,detrend=detrend,boundary='none', &
            padded=.false.,scaling=scale,specmode=.true., &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            spectral_amplitude=chosen_mode /= 'psd', detrend_callback=detrend_callback)
        r%frequencies = temp%frequencies
        r%times = temp%times
        if (temp%status == spectral_status_success) then
            if (chosen_mode == 'complex') then
                r%values = temp%values
                allocate(r%power(0,0))
            else
                select case (chosen_mode)
                case ('psd')
                    r%power = real(temp%values,dp)
                case ('magnitude')
                    r%power = abs(temp%values)
                case ('angle', 'phase')
                    r%power = atan2(aimag(temp%values),real(temp%values,dp))
                    if (chosen_mode == 'phase') call spectral_unwrap_phase(r%power)
                end select
            end if
        else
            allocate(r%power(0,0))
        end if
        r%status = temp%status
    end function spectrogram_complex

    function compute_complex_stft(x,fs,window,nperseg,noverlap,nfft,detrend, &
            boundary,padded,scaling,specmode, &
            window_parameter, window_parameter2, window_coefficients, window_values, spectral_amplitude, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! unpadded complex signal
        real(dp), intent(in), optional :: fs !! sampling rate
        character(*), intent(in), optional :: window !! analysis window
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap size
        integer, intent(in), optional :: nfft !! DFT length
        character(*), intent(in), optional :: detrend !! detrend mode
        character(*), intent(in), optional :: boundary !! zeros, even, odd, constant, or none
        logical, intent(in), optional :: padded !! pad the final frame
        character(*), intent(in), optional :: scaling !! scale choice
        logical, intent(in) :: specmode !! squared magnitude for spectrogram

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        logical, intent(in), optional :: spectral_amplitude !! complex STFT instead of PSD for spectrogram modes
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(stft_result) :: r
        real(dp), allocatable :: w(:), effective_parameter
        complex(dp), allocatable :: signal(:), ft(:)
        real(dp) :: rate, factor, powerfactor
        integer :: nw,no,nf,left,nall,hop,nframes,total,j
        logical :: ok,padmore,want_amplitude
        character(:), allocatable :: win, trend, edges, scale

        want_amplitude = .false.
        if (present(spectral_amplitude)) want_amplitude = spectral_amplitude
        allocate(r%frequencies(0),r%times(0),r%values(0,0))
        if (size(x) < 1) then
            r%status = spectral_status_success
            return
        end if
        if (.not. finite_complex(x)) return
        win = 'hann'
        if (specmode) win = 'tukey'
        if (present(window)) win = trim(window)
        if (present(window_parameter)) effective_parameter = window_parameter
        if (specmode .and. .not. present(window) .and. .not. present(window_parameter) .and. &
            .not. present(window_values) .and. .not. present(window_coefficients)) &
            effective_parameter = 0.25_dp
        scale = 'spectrum'
        if (specmode) scale = 'density'
        if (present(scaling)) scale = trim(scaling)
        if (specmode) then
            call spectral_setup(size(x),fs,window=window,nperseg=nperseg, &
                noverlap=noverlap,nfft=nfft,scaling=scale,default_window=win, &
                default_overlap=8,nw=nw,no=no,nf=nf,rate=rate,w=w, &
                factor=powerfactor,ok=ok, &
                window_parameter=effective_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients, window_values=window_values)
            if (.not. ok) return
            factor = sqrt(powerfactor)
        else
            if (scale == 'spectrum') then
                call spectral_setup(size(x),fs,window=window,nperseg=nperseg, &
                    noverlap=noverlap,nfft=nfft,scaling=scale,default_window=win, &
                    default_overlap=0,nw=nw,no=no,nf=nf,rate=rate,w=w, &
                    factor=powerfactor,ok=ok, &
                window_parameter=effective_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients, window_values=window_values)
            else if (scale == 'psd') then
                call spectral_setup(size(x),fs,window=window,nperseg=nperseg, &
                    noverlap=noverlap,nfft=nfft,scaling='density',default_window=win, &
                    default_overlap=0,nw=nw,no=no,nf=nf,rate=rate,w=w, &
                    factor=powerfactor,ok=ok, &
                window_parameter=effective_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients, window_values=window_values)
            else
                return
            end if
            if (.not. ok) return
            factor = sqrt(powerfactor)
        end if
        trend = 'none'
        if (specmode) trend = 'constant'
        if (present(detrend)) trend = trim(detrend)
        if (present(detrend_callback) .and. present(detrend)) return
        if (.not. present(detrend_callback)) then
            if (trend /= 'none' .and. trend /= 'constant' .and. &
                trend /= 'linear' .and. trend /= 'false') return
        end if
        edges = 'zeros'
        if (specmode) edges = 'none'
        if (present(boundary)) edges = trim(boundary)
        if (.not. stft_edge_mode_valid(edges)) return
        padmore = .true.
        if (specmode) padmore = .false.
        if (present(padded)) padmore = padded
        left = 0
        if (edges /= 'none') left = nw/2
        nall = size(x)+2*left
        hop = nw-no
        if (padmore) then
            nframes = 1+(nall-nw+hop-1)/hop
        else
            nframes = 1+(nall-nw)/hop
        end if
        nframes = max(1,nframes)
        total = nw+(nframes-1)*hop
        allocate(signal(max(total,nall)),ft(nf))
        signal = (0.0_dp,0.0_dp)
        call fill_stft_edges(x,left,edges,signal)
        deallocate(r%frequencies,r%times,r%values)
        allocate(r%frequencies(nf),r%times(nframes),r%values(nf,nframes))
        call two_sided_frequencies(nf,rate,r%frequencies)
        do j = 1, nframes
            r%times(j) = (real((j-1)*hop,dp)+real(nw,dp)/2.0_dp - &
                merge(real(nw,dp)/2.0_dp,0.0_dp,edges/='none'))/rate
            call transform_complex_segment(signal,1+(j-1)*hop,nw,nf,w,trend,ft,ok,detrend_callback)
            if (.not. ok) return
            if (specmode .and. .not. want_amplitude) then
                r%values(:,j) = cmplx(abs(ft)**2*powerfactor,0.0_dp,dp)
            else
                r%values(:,j) = ft*factor
            end if
        end do
        r%status = spectral_status_success
    end function compute_complex_stft


    ! Mixed real/complex spectra use a two-sided FFT, matching the spectrum
    ! obtained when the real operand is explicitly promoted to complex(dp).
    ! The legacy SciPy 1.17 dispatcher has order-dependent mixed-dtype errors;
    ! the explicitly promoted calculation is used as the numerical reference.

    function csd_real_complex(x, y, fs, window, nperseg, noverlap, nfft, detrend, &
            scaling, average, return_onesided, window_parameter, window_parameter2, &
            window_coefficients, window_values, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! first input signal
        complex(dp), intent(in) :: y(:) !! second input signal
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! analysis window name
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! shared samples between segments
        integer, intent(in), optional :: nfft !! transform length
        character(*), intent(in), optional :: detrend !! constant, linear, none
        character(*), intent(in), optional :: scaling !! density or spectrum
        character(*), intent(in), optional :: average !! mean or median
        logical, intent(in), optional :: return_onesided !! ignored when either input is complex
        real(dp), intent(in), optional :: window_parameter !! first window parameter
        real(dp), intent(in), optional :: window_parameter2 !! second window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit periodic window
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! complex segment callback
        type(cross_spectral_result) :: r

        r = csd_complex(cmplx(x,0.0_dp,dp), y, fs=fs, window=window, &
            nperseg=nperseg, noverlap=noverlap, nfft=nfft, detrend=detrend, &
            scaling=scaling, average=average, return_onesided=return_onesided, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
    end function csd_real_complex

    function csd_complex_real(x, y, fs, window, nperseg, noverlap, nfft, detrend, &
            scaling, average, return_onesided, window_parameter, window_parameter2, &
            window_coefficients, window_values, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! first input signal
        real(dp), intent(in) :: y(:) !! second input signal
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! analysis window name
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! shared samples between segments
        integer, intent(in), optional :: nfft !! transform length
        character(*), intent(in), optional :: detrend !! constant, linear, none
        character(*), intent(in), optional :: scaling !! density or spectrum
        character(*), intent(in), optional :: average !! mean or median
        logical, intent(in), optional :: return_onesided !! ignored when either input is complex
        real(dp), intent(in), optional :: window_parameter !! first window parameter
        real(dp), intent(in), optional :: window_parameter2 !! second window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general-cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit periodic window
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! complex segment callback
        type(cross_spectral_result) :: r

        r = csd_complex(x, cmplx(y,0.0_dp,dp), fs=fs, window=window, &
            nperseg=nperseg, noverlap=noverlap, nfft=nfft, detrend=detrend, &
            scaling=scaling, average=average, return_onesided=return_onesided, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
    end function csd_complex_real

    function coherence_real_complex(x, y, fs, window, nperseg, noverlap, nfft, detrend, &
            window_parameter, window_parameter2, window_coefficients, &
            window_values, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! first input signal
        complex(dp), intent(in) :: y(:) !! second input signal
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! analysis window name
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlapping samples
        integer, intent(in), optional :: nfft !! transform length
        character(*), intent(in), optional :: detrend !! constant, linear, none
        real(dp), intent(in), optional :: window_parameter !! first window parameter
        real(dp), intent(in), optional :: window_parameter2 !! second window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine weights
        real(dp), intent(in), optional :: window_values(:) !! explicit window samples
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! complex segment callback
        type(spectral_result) :: r

        r = coherence_complex(cmplx(x,0.0_dp,dp), y, fs=fs, window=window, &
            nperseg=nperseg, noverlap=noverlap, nfft=nfft, detrend=detrend, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
    end function coherence_real_complex

    function coherence_complex_real(x, y, fs, window, nperseg, noverlap, nfft, detrend, &
            window_parameter, window_parameter2, window_coefficients, &
            window_values, detrend_callback) result(r)
        complex(dp), intent(in) :: x(:) !! first input signal
        real(dp), intent(in) :: y(:) !! second input signal
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! analysis window name
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlapping samples
        integer, intent(in), optional :: nfft !! transform length
        character(*), intent(in), optional :: detrend !! constant, linear, none
        real(dp), intent(in), optional :: window_parameter !! first window parameter
        real(dp), intent(in), optional :: window_parameter2 !! second window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine weights
        real(dp), intent(in), optional :: window_values(:) !! explicit window samples
        procedure(spectral_complex_detrend_callback), optional :: detrend_callback !! complex segment callback
        type(spectral_result) :: r

        r = coherence_complex(x, cmplx(y,0.0_dp,dp), fs=fs, window=window, &
            nperseg=nperseg, noverlap=noverlap, nfft=nfft, detrend=detrend, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
    end function coherence_complex_real

end module scifort_signal_spectral_complex
