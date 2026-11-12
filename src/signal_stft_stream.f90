! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Streaming forward STFT and power spectrograms with bounded persistent storage.
! Explicit nperseg is required because a stream's final length is not known.
module scifort_signal_stft_stream
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, forward_fft_planned, get_window, &
        spectral_status_success, spectral_status_invalid_input
    implicit none
    private

    type, public :: stft_stream
        private
        integer :: nperseg = 0, nfft = 0, nfreq = 0, hop = 0, left = 0
        integer(int64) :: total = 0_int64, received = 0_int64, next_start = 0_int64
        real(dp) :: fs = 1.0_dp, factor = 1.0_dp
        logical :: active = .false., finished = .false., input_complex = .false.
        logical :: onesided = .true., padded = .true., power_mode = .false.
        logical :: boundary_zeros = .true.
        character(8) :: detrend = 'none'
        real(dp), allocatable :: window(:)
        complex(dp), allocatable :: ring(:), work(:)
        type(fft_plan) :: transform_plan
    end type stft_stream

    interface stft_stream_push
        module procedure stft_push_real, stft_push_complex
    end interface
    interface spectrogram_stream_push
        module procedure spectrogram_push_real, spectrogram_push_complex
    end interface
    public :: stft_stream_init, stft_stream_push, stft_stream_finish
    public :: spectrogram_stream_init, spectrogram_stream_push, spectrogram_stream_finish
    public :: stft_stream_frequencies

contains

    subroutine stft_stream_init(state, nperseg, status, fs, window, noverlap, nfft, &
        detrend, boundary, padded, scaling, window_values, input_complex, return_onesided)
        type(stft_stream), intent(inout) :: state !! state, unchanged if initialization fails
        integer, intent(in) :: nperseg !! explicit analysis-window length
        integer, intent(out) :: status !! spectral status code
        real(dp), intent(in), optional :: fs !! samples per second
        character(*), intent(in), optional :: window !! default 'hann'
        integer, intent(in), optional :: noverlap !! default half a segment
        integer, intent(in), optional :: nfft !! FFT length, at least nperseg
        character(*), intent(in), optional :: detrend !! none, constant, or linear
        character(*), intent(in), optional :: boundary !! 'zeros' or 'none'
        logical, intent(in), optional :: padded !! pad final segment, default true
        character(*), intent(in), optional :: scaling !! 'spectrum' or 'psd'
        real(dp), intent(in), optional :: window_values(:) !! custom periodic window values
        logical, intent(in), optional :: input_complex !! true for complex input samples
        logical, intent(in), optional :: return_onesided !! default true for real input

        call initialize(state,nperseg,status,.false.,fs,window,noverlap,nfft, &
            detrend,boundary,padded,scaling,window_values,input_complex,return_onesided)
    end subroutine stft_stream_init

    subroutine spectrogram_stream_init(state, nperseg, status, fs, window, noverlap, &
        nfft, detrend, scaling, window_values, input_complex, return_onesided)
        type(stft_stream), intent(inout) :: state !! state, unchanged if initialization fails
        integer, intent(in) :: nperseg !! explicit analysis-window length
        integer, intent(out) :: status !! spectral status code
        real(dp), intent(in), optional :: fs !! samples per second
        character(*), intent(in), optional :: window !! default periodic Tukey with alpha 0.25
        integer, intent(in), optional :: noverlap !! default nperseg/8
        integer, intent(in), optional :: nfft !! transform length
        character(*), intent(in), optional :: detrend !! default 'constant'
        character(*), intent(in), optional :: scaling !! 'density' or 'spectrum'
        real(dp), intent(in), optional :: window_values(:) !! user-supplied window
        logical, intent(in), optional :: input_complex !! true for complex input samples
        logical, intent(in), optional :: return_onesided !! default true for real input

        call initialize(state,nperseg,status,.true.,fs,window,noverlap,nfft, &
            detrend,scaling=scaling,window_values=window_values, &
            input_complex=input_complex,return_onesided=return_onesided)
    end subroutine spectrogram_stream_init

    subroutine initialize(state,nperseg,status,power_mode,fs,window,noverlap,nfft, &
        detrend,boundary,padded,scaling,window_values,input_complex,return_onesided)
        implicit none
        type(stft_stream), intent(inout) :: state !! stream being prepared
        integer, intent(in) :: nperseg !! positive segment length
        integer, intent(out) :: status !! spectral status code
        logical, intent(in) :: power_mode !! whether output is a power spectrogram
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! periodic window name
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: nfft !! transform length
        character(*), intent(in), optional :: detrend !! detrending mode
        character(*), intent(in), optional :: boundary !! zeros or none
        logical, intent(in), optional :: padded !! pad final segment
        character(*), intent(in), optional :: scaling !! normalization choice
        real(dp), intent(in), optional :: window_values(:) !! explicit window weights
        logical, intent(in), optional :: input_complex !! complex input flag
        logical, intent(in), optional :: return_onesided !! one-sided real transform
        type(stft_stream) :: trial
        integer :: no, nf, ierr, j
        real(dp) :: energy, norm
        character(:), allocatable :: name, trend, edges, scale

        status = spectral_status_invalid_input
        if (nperseg < 1) return
        if (nperseg > huge(nperseg)/4) return
        trial%nperseg = nperseg
        no = nperseg/2
        if (power_mode) no = nperseg/8
        if (present(noverlap)) no = noverlap
        if (no < 0 .or. no >= nperseg) return
        trial%hop = nperseg-no
        nf = nperseg
        if (present(nfft)) nf = nfft
        if (nf < nperseg) return
        trial%nfft = nf
        trial%fs = 1.0_dp
        if (present(fs)) trial%fs = fs
        if (.not. ieee_is_finite(trial%fs)) return
        if (trial%fs <= 0.0_dp) return
        trial%input_complex = .false.
        if (present(input_complex)) trial%input_complex = input_complex
        trial%onesided = .not.trial%input_complex
        if (present(return_onesided)) trial%onesided = return_onesided
        if (trial%input_complex) trial%onesided = .false.
        trial%nfreq = nf
        if (trial%onesided) trial%nfreq = nf/2+1
        trial%power_mode = power_mode
        trial%padded = .not.power_mode
        if (present(padded)) trial%padded = padded
        edges = 'zeros'
        if (power_mode) edges = 'none'
        if (present(boundary)) edges = trim(boundary)
        if (edges /= 'zeros' .and. edges /= 'none') return
        trial%boundary_zeros = edges == 'zeros'
        trial%left = 0
        if (trial%boundary_zeros) trial%left = nperseg/2
        trend = 'none'
        if (power_mode) trend = 'constant'
        if (present(detrend)) trend = trim(detrend)
        if (trend /= 'none' .and. trend /= 'false' .and. &
            trend /= 'constant' .and. trend /= 'linear') return
        trial%detrend = trend
        if (present(window) .and. present(window_values)) return
        if (present(window_values)) then
            if (size(window_values) /= nperseg) return
            trial%window = window_values
        else
            name = 'hann'
            if (power_mode) name = 'tukey'
            if (present(window)) name = trim(window)
            if (power_mode .and. .not. present(window)) then
                ! Match the batch spectrogram default; a named Tukey still uses alpha 0.5.
                trial%window = get_window(name,nperseg,parameter=0.25_dp,status=ierr)
            else
                trial%window = get_window(name,nperseg,status=ierr)
            end if
            if (ierr /= spectral_status_success) return
        end if
        if (.not. all(ieee_is_finite(trial%window))) return
        energy = sum(trial%window*trial%window)
        if (.not. ieee_is_finite(energy)) return
        scale = 'spectrum'
        if (power_mode) scale = 'density'
        if (present(scaling)) scale = trim(scaling)
        if (power_mode) then
            select case(scale)
            case('density')
                norm = trial%fs*energy
            case('spectrum')
                norm = sum(trial%window)**2
            case default
                return
            end select
            if (.not. ieee_is_finite(norm) .or. norm <= tiny(norm)) return
            trial%factor = 1.0_dp/norm
        else
            select case(scale)
            case('spectrum')
                norm = sum(trial%window)
                if (.not. ieee_is_finite(norm) .or. abs(norm) <= sqrt(tiny(norm))) return
                trial%factor = 1.0_dp/norm
            case('psd')
                norm = trial%fs*energy
                if (.not. ieee_is_finite(norm) .or. norm <= tiny(norm)) return
                trial%factor = 1.0_dp/sqrt(norm)
            case default
                return
            end select
        end if
        if (.not. ieee_is_finite(trial%factor)) return
        allocate(trial%ring(nperseg),trial%work(nf))
        trial%ring = (0.0_dp,0.0_dp)
        trial%active = .true.
        ! Logical left-edge zeros are already present, before the first push.
        do j=1,trial%left
            trial%ring(j) = (0.0_dp,0.0_dp)
        end do
        trial%total = int(trial%left,int64)
        call trial%transform_plan%initialize(nf, ierr)
        if (ierr /= spectral_status_success) return
        state = trial
        status = spectral_status_success
    end subroutine initialize

    subroutine stft_stream_frequencies(state, frequencies, status)
        type(stft_stream), intent(in) :: state !! initialized stream
        real(dp), allocatable, intent(out) :: frequencies(:) !! one- or two-sided FFT bins
        integer, intent(out) :: status !! spectral status code
        integer :: k, signed_k

        allocate(frequencies(0))
        status = spectral_status_invalid_input
        if (.not.state%active) return
        deallocate(frequencies)
        allocate(frequencies(state%nfreq))
        do k=0,state%nfreq-1
            signed_k = k
            if (.not.state%onesided .and. k > (state%nfft-1)/2) signed_k = k-state%nfft
            frequencies(k+1) = real(signed_k,dp)*state%fs/real(state%nfft,dp)
        end do
        status = spectral_status_success
    end subroutine stft_stream_frequencies

    subroutine transform_frame(state,spectrum)
        type(stft_stream), intent(inout) :: state !! initialized state
        complex(dp), intent(out) :: spectrum(:) !! current frequency-domain frame
        integer :: i, j
        complex(dp) :: avg, slope
        real(dp) :: midpoint, denominator

        state%work = (0.0_dp,0.0_dp)
        do i=1,state%nperseg
            j = int(mod(state%next_start+int(i-1,int64),int(state%nperseg,int64)))+1
            state%work(i) = state%ring(j)
        end do
        select case(state%detrend)
        case('constant')
            avg = sum(state%work(1:state%nperseg))/real(state%nperseg,dp)
            state%work(1:state%nperseg) = state%work(1:state%nperseg)-avg
        case('linear')
            avg = sum(state%work(1:state%nperseg))/real(state%nperseg,dp)
            midpoint = (real(state%nperseg,dp)+1.0_dp)*0.5_dp
            denominator = real(state%nperseg,dp)* &
                (real(state%nperseg,dp)**2-1.0_dp)/12.0_dp
            slope = (0.0_dp,0.0_dp)
            if (state%nperseg > 1) then
                do i=1,state%nperseg
                    slope = slope + (real(i,dp)-midpoint)*(state%work(i)-avg)
                end do
                slope = slope/denominator
            end if
            do i=1,state%nperseg
                state%work(i) = state%work(i)-avg-slope*(real(i,dp)-midpoint)
            end do
        end select
        state%work(1:state%nperseg) = state%work(1:state%nperseg)*state%window
        call forward_fft_planned(state%work, state%transform_plan)
        if (state%power_mode) then
            spectrum = cmplx(abs(state%work(1:state%nfreq))**2*state%factor,0.0_dp,dp)
            if (state%onesided) then
                do i=2,state%nfreq
                    if (mod(state%nfft,2)==0 .and. i==state%nfreq) cycle
                    spectrum(i) = 2.0_dp*spectrum(i)
                end do
            end if
        else
            spectrum = state%work(1:state%nfreq)*state%factor
        end if
    end subroutine transform_frame

    subroutine accept_sample(state,x,spectra,times,k)
        type(stft_stream), intent(inout) :: state !! stream being advanced
        complex(dp), intent(in) :: x !! next time-domain sample
        complex(dp), intent(inout) :: spectra(:,:) !! frame output space
        real(dp), intent(inout) :: times(:) !! frame time-coordinate space
        integer, intent(inout) :: k !! emitted frame count for this call
        integer :: index

        index = int(mod(state%total,int(state%nperseg,int64)))+1
        state%ring(index) = x
        state%total = state%total+1_int64
        if (state%total >= state%next_start+int(state%nperseg,int64)) then
            k = k+1
            call transform_frame(state,spectra(:,k))
            times(k) = real(state%next_start,dp)/state%fs
            if (.not.state%boundary_zeros) then
                times(k) = (real(state%next_start,dp)+real(state%nperseg,dp)/2.0_dp)/state%fs
            end if
            state%next_start = state%next_start+int(state%hop,int64)
        end if
    end subroutine accept_sample

    subroutine prepare_push(state,n,complex_input,mode,spectra,times,status)
        type(stft_stream), intent(in) :: state !! stream state, not modified
        integer, intent(in) :: n !! samples in incoming chunk
        logical, intent(in) :: complex_input !! type of incoming samples
        logical, intent(in) :: mode !! true for power output
        complex(dp), allocatable, intent(out) :: spectra(:,:) !! allocated output buffer
        real(dp), allocatable, intent(out) :: times(:) !! allocated time buffer
        integer, intent(out) :: status !! spectral status code
        integer(int64) :: upto, frames

        status = spectral_status_invalid_input
        allocate(spectra(0,0),times(0))
        if (.not.state%active .or. state%finished) return
        if (state%input_complex .neqv. complex_input) return
        if (state%power_mode .neqv. mode) return
        if (state%total > huge(upto)-int(n,int64)-2_int64*int(state%nperseg,int64)) return
        upto = state%total+int(n,int64)
        frames = 0_int64
        if (upto >= state%next_start+int(state%nperseg,int64)) then
            frames = 1_int64+(upto-state%next_start-int(state%nperseg,int64))/int(state%hop,int64)
        end if
        if (frames > int(huge(n),int64)) return
        deallocate(spectra,times)
        allocate(spectra(state%nfreq,int(frames)),times(int(frames)))
        status = spectral_status_success
    end subroutine prepare_push

    subroutine push_real_core(state,x,spectra,times,status,mode)
        type(stft_stream), intent(inout) :: state !! stream advanced only on success
        real(dp), intent(in) :: x(:) !! time-domain data chunk
        complex(dp), allocatable, intent(out) :: spectra(:,:) !! new STFT frames
        real(dp), allocatable, intent(out) :: times(:) !! new frame times
        integer, intent(out) :: status !! success or rejected chunk
        logical, intent(in) :: mode !! whether this is a spectrogram call
        integer :: i, k

        call prepare_push(state,size(x),.false.,mode,spectra,times,status)
        if (status /= spectral_status_success) return
        if (.not.all(ieee_is_finite(x))) then
            status = spectral_status_invalid_input
            deallocate(spectra,times)
            allocate(spectra(0,0),times(0))
            return
        end if
        k = 0
        do i=1,size(x)
            call accept_sample(state,cmplx(x(i),0.0_dp,dp),spectra,times,k)
        end do
        state%received = state%received+int(size(x),int64)
    end subroutine push_real_core

    subroutine push_complex_core(state,x,spectra,times,status,mode)
        type(stft_stream), intent(inout) :: state !! stream advanced only on success
        complex(dp), intent(in) :: x(:) !! complex data chunk
        complex(dp), allocatable, intent(out) :: spectra(:,:) !! new STFT frames
        real(dp), allocatable, intent(out) :: times(:) !! new frame times
        integer, intent(out) :: status !! success or rejected chunk
        logical, intent(in) :: mode !! whether this is a spectrogram call
        integer :: i, k

        call prepare_push(state,size(x),.true.,mode,spectra,times,status)
        if (status /= spectral_status_success) return
        if (.not.all(ieee_is_finite(real(x,dp))) .or. &
            .not.all(ieee_is_finite(aimag(x)))) then
            status = spectral_status_invalid_input
            deallocate(spectra,times)
            allocate(spectra(0,0),times(0))
            return
        end if
        k = 0
        do i=1,size(x)
            call accept_sample(state,x(i),spectra,times,k)
        end do
        state%received = state%received+int(size(x),int64)
    end subroutine push_complex_core

    subroutine stft_push_real(state,x,values,status,times)
        type(stft_stream), intent(inout) :: state !! initialized real-input stream
        real(dp), intent(in) :: x(:) !! input chunk
        complex(dp), allocatable, intent(out) :: values(:,:) !! new STFT frames
        integer, intent(out) :: status !! success or invalid input
        real(dp), allocatable, intent(out), optional :: times(:) !! frame times
        real(dp), allocatable :: chunk_times(:)

        call push_real_core(state,x,values,chunk_times,status,.false.)
        if (present(times)) times = chunk_times
    end subroutine stft_push_real

    subroutine stft_push_complex(state,x,values,status,times)
        type(stft_stream), intent(inout) :: state !! initialized complex-input stream
        complex(dp), intent(in) :: x(:) !! input chunk
        complex(dp), allocatable, intent(out) :: values(:,:) !! new STFT frames
        integer, intent(out) :: status !! success or invalid input
        real(dp), allocatable, intent(out), optional :: times(:) !! frame times
        real(dp), allocatable :: chunk_times(:)

        call push_complex_core(state,x,values,chunk_times,status,.false.)
        if (present(times)) times = chunk_times
    end subroutine stft_push_complex

    subroutine finish_core(state,values,times,status,power_mode)
        type(stft_stream), intent(inout) :: state !! stream to finalize
        complex(dp), allocatable, intent(out) :: values(:,:) !! final frames
        real(dp), allocatable, intent(out) :: times(:) !! final frame times
        integer, intent(out) :: status !! spectral status code
        logical, intent(in) :: power_mode !! desired stream mode
        integer(int64) :: nall, final_size, nframes, zeros, j
        integer :: k

        allocate(values(0,0),times(0))
        status = spectral_status_invalid_input
        if (.not.state%active .or. state%finished) return
        if (state%power_mode .neqv. power_mode) return
        if (state%received == 0_int64) then
            state%finished = .true.
            status = spectral_status_success
            return
        end if
        nall = state%received+2_int64*int(state%left,int64)
        if (state%padded) then
            nframes = max(1_int64,1_int64+ &
                (max(0_int64,nall-int(state%nperseg,int64))+int(state%hop,int64)-1_int64)/int(state%hop,int64))
        else
            nframes = max(1_int64,1_int64+ &
                (nall-int(state%nperseg,int64))/int(state%hop,int64))
        end if
        final_size = max(nall,int(state%nperseg,int64)+(nframes-1_int64)*int(state%hop,int64))
        zeros = final_size-state%total
        if (zeros < 0_int64) return
        nframes = 0_int64
        if (final_size >= state%next_start+int(state%nperseg,int64)) then
            nframes = 1_int64+(final_size-state%next_start-int(state%nperseg,int64))/int(state%hop,int64)
        end if
        if (nframes > int(huge(k),int64)) return
        deallocate(values,times)
        allocate(values(state%nfreq,int(nframes)),times(int(nframes)))
        k = 0
        do j=1_int64,zeros
            call accept_sample(state,cmplx(0.0_dp,0.0_dp,dp),values,times,k)
        end do
        state%finished = .true.
        status = spectral_status_success
    end subroutine finish_core

    subroutine stft_stream_finish(state,values,status,times)
        type(stft_stream), intent(inout) :: state !! active analysis stream
        complex(dp), allocatable, intent(out) :: values(:,:) !! final frames
        integer, intent(out) :: status !! success or invalid input
        real(dp), allocatable, intent(out), optional :: times(:) !! final frame times
        real(dp), allocatable :: chunk_times(:)

        call finish_core(state,values,chunk_times,status,.false.)
        if (present(times)) times = chunk_times
    end subroutine stft_stream_finish

    subroutine spectrogram_push_real(state,x,power,status,times)
        type(stft_stream), intent(inout) :: state !! active real spectrogram stream
        real(dp), intent(in) :: x(:) !! real input chunk
        real(dp), allocatable, intent(out) :: power(:,:) !! new PSD frames
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, intent(out), optional :: times(:) !! frame times
        complex(dp), allocatable :: values(:,:)
        real(dp), allocatable :: chunk_times(:)

        call push_real_core(state,x,values,chunk_times,status,.true.)
        power = real(values,dp)
        if (present(times)) times = chunk_times
    end subroutine spectrogram_push_real

    subroutine spectrogram_push_complex(state,x,power,status,times)
        type(stft_stream), intent(inout) :: state !! active complex spectrogram stream
        complex(dp), intent(in) :: x(:) !! complex input chunk
        real(dp), allocatable, intent(out) :: power(:,:) !! new PSD frames
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, intent(out), optional :: times(:) !! frame times
        complex(dp), allocatable :: values(:,:)
        real(dp), allocatable :: chunk_times(:)

        call push_complex_core(state,x,values,chunk_times,status,.true.)
        power = real(values,dp)
        if (present(times)) times = chunk_times
    end subroutine spectrogram_push_complex

    subroutine spectrogram_stream_finish(state,power,status,times)
        type(stft_stream), intent(inout) :: state !! active spectrogram stream
        real(dp), allocatable, intent(out) :: power(:,:) !! final PSD frames
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, intent(out), optional :: times(:) !! final frame times
        complex(dp), allocatable :: values(:,:)
        real(dp), allocatable :: chunk_times(:)

        call finish_core(state,values,chunk_times,status,.true.)
        power = real(values,dp)
        if (present(times)) times = chunk_times
    end subroutine spectrogram_stream_finish

end module scifort_signal_stft_stream
