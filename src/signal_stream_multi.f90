! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Synchronous multichannel streaming STFT, spectrogram and inverse STFT.
! Inputs are (sample, channel); forward outputs are (frequency, channel, frame).
! Inverse inputs are (frequency, channel), with outputs (sample, channel).
! Only fixed channel counts are supported within a stream.
module scifort_signal_stream_multi
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success, spectral_status_invalid_input
    use scifort_signal_istft, only : istft_status_success, istft_status_invalid_input, &
        istft_status_noninvertible
    use scifort_signal_stft_stream, only : stft_stream, stft_stream_init, &
        stft_stream_push, stft_stream_finish, stft_stream_frequencies, &
        spectrogram_stream_init, spectrogram_stream_push, spectrogram_stream_finish
    use scifort_signal_istft_stream, only : istft_stream, istft_stream_init, &
        istft_stream_push, istft_stream_finish
    implicit none
    private

    type, public :: stft_stream_multi
        private
        type(stft_stream), allocatable :: channels(:)
        integer :: nchannels = 0
        logical :: active = .false., finished = .false.
        logical :: power_mode = .false., input_complex = .false.
    end type stft_stream_multi

    type, public :: istft_stream_multi
        private
        type(istft_stream), allocatable :: channels(:)
        integer :: nchannels = 0, nfreq = 0
        logical :: active = .false., finished = .false., onesided = .true.
    end type istft_stream_multi

    interface stft_stream_multi_push
        module procedure forward_push_real, forward_push_complex
    end interface
    interface spectrogram_stream_multi_push
        module procedure spectrogram_push_real, spectrogram_push_complex
    end interface
    interface istft_stream_multi_push
        module procedure inverse_push_real, inverse_push_complex
    end interface
    interface istft_stream_multi_finish
        module procedure inverse_finish_real, inverse_finish_complex
    end interface

    public :: stft_stream_multi_init, stft_stream_multi_push, stft_stream_multi_finish
    public :: stft_stream_multi_frequencies
    public :: spectrogram_stream_multi_init, spectrogram_stream_multi_push, &
        spectrogram_stream_multi_finish
    public :: istft_stream_multi_init, istft_stream_multi_push, istft_stream_multi_finish

contains

    subroutine stft_stream_multi_init(state, nchannels, nperseg, status, fs, window, &
        noverlap, nfft, detrend, boundary, padded, scaling, window_values, &
        input_complex, return_onesided)
        type(stft_stream_multi), intent(inout) :: state !! unchanged on invalid configuration
        integer, intent(in) :: nchannels !! number of independent signals, positive
        integer, intent(in) :: nperseg !! analysis window length
        integer, intent(out) :: status !! spectral status code
        real(dp), optional, intent(in) :: fs !! sampling rate
        character(*), optional, intent(in) :: window !! named analysis window
        integer, optional, intent(in) :: noverlap !! overlap samples
        integer, optional, intent(in) :: nfft !! transform length
        character(*), optional, intent(in) :: detrend !! detrending mode
        character(*), optional, intent(in) :: boundary !! 'zeros' or 'none'
        logical, optional, intent(in) :: padded !! pad final frame
        character(*), optional, intent(in) :: scaling !! 'spectrum' or 'psd'
        real(dp), optional, intent(in) :: window_values(:) !! explicit window
        logical, optional, intent(in) :: input_complex !! complex chunk mode
        logical, optional, intent(in) :: return_onesided !! one-sided spectrum for real inputs
        type(stft_stream_multi) :: trial
        integer :: j

        status = spectral_status_invalid_input
        if (nchannels < 1) return
        allocate(trial%channels(nchannels))
        do j=1,nchannels
            call stft_stream_init(trial%channels(j),nperseg,status,fs,window, &
                noverlap,nfft,detrend,boundary,padded,scaling,window_values, &
                input_complex,return_onesided)
            if (status /= spectral_status_success) return
        end do
        trial%nchannels = nchannels
        trial%input_complex = .false.
        if (present(input_complex)) trial%input_complex = input_complex
        trial%active = .true.
        state = trial
    end subroutine stft_stream_multi_init

    subroutine spectrogram_stream_multi_init(state, nchannels, nperseg, status, fs, &
        window, noverlap, nfft, detrend, scaling, window_values, input_complex, &
        return_onesided)
        type(stft_stream_multi), intent(inout) :: state !! unchanged on invalid configuration
        integer, intent(in) :: nchannels !! number of independent signals, positive
        integer, intent(in) :: nperseg !! analysis window length
        integer, intent(out) :: status !! spectral status code
        real(dp), optional, intent(in) :: fs !! sampling rate
        character(*), optional, intent(in) :: window !! named analysis window
        integer, optional, intent(in) :: noverlap !! overlap samples
        integer, optional, intent(in) :: nfft !! transform length
        character(*), optional, intent(in) :: detrend !! detrending mode
        character(*), optional, intent(in) :: scaling !! 'density' or 'spectrum'
        real(dp), optional, intent(in) :: window_values(:) !! explicit window
        logical, optional, intent(in) :: input_complex !! complex chunk mode
        logical, optional, intent(in) :: return_onesided !! one-sided spectrum for real inputs
        type(stft_stream_multi) :: trial
        integer :: j

        status = spectral_status_invalid_input
        if (nchannels < 1) return
        allocate(trial%channels(nchannels))
        do j=1,nchannels
            call spectrogram_stream_init(trial%channels(j),nperseg,status,fs,window, &
                noverlap,nfft,detrend,scaling,window_values,input_complex,return_onesided)
            if (status /= spectral_status_success) return
        end do
        trial%nchannels = nchannels
        trial%input_complex = .false.
        if (present(input_complex)) trial%input_complex = input_complex
        trial%power_mode = .true.
        trial%active = .true.
        state = trial
    end subroutine spectrogram_stream_multi_init

    subroutine stft_stream_multi_frequencies(state, frequencies, status)
        type(stft_stream_multi), intent(in) :: state !! initialized multichannel stream
        real(dp), allocatable, intent(out) :: frequencies(:) !! common frequency bins
        integer, intent(out) :: status !! success or invalid input

        allocate(frequencies(0))
        status = spectral_status_invalid_input
        if (.not.state%active) return
        call stft_stream_frequencies(state%channels(1),frequencies,status)
    end subroutine stft_stream_multi_frequencies

    subroutine forward_real_core(state,x,values,times,status,power_mode)
        type(stft_stream_multi), intent(inout) :: state !! active multichannel stream
        real(dp), intent(in) :: x(:,:) !! incoming (sample,channel) chunk
        complex(dp), allocatable, intent(out) :: values(:,:,:) !! (frequency,channel,frame)
        real(dp), allocatable, intent(out) :: times(:) !! common frame coordinates
        integer, intent(out) :: status !! spectral status code
        logical, intent(in) :: power_mode !! return power instead of STFT
        complex(dp), allocatable :: one(:,:)
        real(dp), allocatable :: power(:,:), t(:)
        integer :: j, nf, nt

        allocate(values(0,0,0),times(0))
        status = spectral_status_invalid_input
        if (.not.state%active .or. state%finished) return
        if (state%power_mode .neqv. power_mode) return
        if (state%input_complex) return
        if (size(x,2) /= state%nchannels) return
        if (.not.all(ieee_is_finite(x))) return
        do j=1,state%nchannels
            if (power_mode) then
                call spectrogram_stream_push(state%channels(j),x(:,j),power,status,t)
                one = cmplx(power,0.0_dp,dp)
            else
                call stft_stream_push(state%channels(j),x(:,j),one,status,t)
            end if
            if (status /= spectral_status_success) return
            if (j==1) then
                nf = size(one,1)
                nt = size(one,2)
                deallocate(values)
                allocate(values(nf,state%nchannels,nt))
                times = t
            else
                if (size(one,1) /= nf .or. size(one,2) /= nt) then
                    status = spectral_status_invalid_input
                    return
                end if
            end if
            values(:,j,:) = one
        end do
    end subroutine forward_real_core

    subroutine forward_complex_core(state,x,values,times,status,power_mode)
        type(stft_stream_multi), intent(inout) :: state !! active multichannel stream
        complex(dp), intent(in) :: x(:,:) !! incoming (sample,channel) chunk
        complex(dp), allocatable, intent(out) :: values(:,:,:) !! (frequency,channel,frame)
        real(dp), allocatable, intent(out) :: times(:) !! common frame coordinates
        integer, intent(out) :: status !! spectral status code
        logical, intent(in) :: power_mode !! return power instead of STFT
        complex(dp), allocatable :: one(:,:)
        real(dp), allocatable :: power(:,:), t(:)
        integer :: j, nf, nt

        allocate(values(0,0,0),times(0))
        status = spectral_status_invalid_input
        if (.not.state%active .or. state%finished) return
        if (state%power_mode .neqv. power_mode) return
        if (.not.state%input_complex) return
        if (size(x,2) /= state%nchannels) return
        if (.not.all(ieee_is_finite(real(x,dp)))) return
        if (.not.all(ieee_is_finite(aimag(x)))) return
        do j=1,state%nchannels
            if (power_mode) then
                call spectrogram_stream_push(state%channels(j),x(:,j),power,status,t)
                one = cmplx(power,0.0_dp,dp)
            else
                call stft_stream_push(state%channels(j),x(:,j),one,status,t)
            end if
            if (status /= spectral_status_success) return
            if (j==1) then
                nf = size(one,1)
                nt = size(one,2)
                deallocate(values)
                allocate(values(nf,state%nchannels,nt))
                times = t
            else
                if (size(one,1) /= nf .or. size(one,2) /= nt) then
                    status = spectral_status_invalid_input
                    return
                end if
            end if
            values(:,j,:) = one
        end do
    end subroutine forward_complex_core

    subroutine forward_push_real(state,x,values,status,times)
        type(stft_stream_multi), intent(inout) :: state !! active multichannel stream
        real(dp), intent(in) :: x(:,:) !! samples by channel
        complex(dp), allocatable, intent(out) :: values(:,:,:) !! frequency by channel by frame
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! shared frame coordinates
        real(dp), allocatable :: t(:)

        call forward_real_core(state,x,values,t,status,.false.)
        if (present(times)) times = t
    end subroutine forward_push_real

    subroutine forward_push_complex(state,x,values,status,times)
        type(stft_stream_multi), intent(inout) :: state !! active multichannel stream
        complex(dp), intent(in) :: x(:,:) !! samples by channel
        complex(dp), allocatable, intent(out) :: values(:,:,:) !! frequency by channel by frame
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! shared frame coordinates
        real(dp), allocatable :: t(:)

        call forward_complex_core(state,x,values,t,status,.false.)
        if (present(times)) times = t
    end subroutine forward_push_complex

    subroutine spectrogram_push_real(state,x,power,status,times)
        type(stft_stream_multi), intent(inout) :: state !! active spectrogram stream
        real(dp), intent(in) :: x(:,:) !! samples by channel
        real(dp), allocatable, intent(out) :: power(:,:,:) !! frequency by channel by frame
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! shared frame coordinates
        real(dp), allocatable :: t(:)
        complex(dp), allocatable :: values(:,:,:)

        call forward_real_core(state,x,values,t,status,.true.)
        power = real(values,dp)
        if (present(times)) times = t
    end subroutine spectrogram_push_real

    subroutine spectrogram_push_complex(state,x,power,status,times)
        type(stft_stream_multi), intent(inout) :: state !! active spectrogram stream
        complex(dp), intent(in) :: x(:,:) !! samples by channel
        real(dp), allocatable, intent(out) :: power(:,:,:) !! frequency by channel by frame
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! shared frame coordinates
        real(dp), allocatable :: t(:)
        complex(dp), allocatable :: values(:,:,:)

        call forward_complex_core(state,x,values,t,status,.true.)
        power = real(values,dp)
        if (present(times)) times = t
    end subroutine spectrogram_push_complex

    subroutine forward_finish_core(state,values,times,status,power_mode)
        type(stft_stream_multi), intent(inout) :: state !! active stream
        complex(dp), allocatable, intent(out) :: values(:,:,:) !! final spectral frames
        real(dp), allocatable, intent(out) :: times(:) !! common frame coordinates
        integer, intent(out) :: status !! spectral status code
        logical, intent(in) :: power_mode !! power spectrogram flag
        complex(dp), allocatable :: one(:,:)
        real(dp), allocatable :: power(:,:), t(:)
        integer :: j, nf, nt

        allocate(values(0,0,0),times(0))
        status = spectral_status_invalid_input
        if (.not.state%active .or. state%finished) return
        if (state%power_mode .neqv. power_mode) return
        do j=1,state%nchannels
            if (power_mode) then
                call spectrogram_stream_finish(state%channels(j),power,status,t)
                one = cmplx(power,0.0_dp,dp)
            else
                call stft_stream_finish(state%channels(j),one,status,t)
            end if
            if (status /= spectral_status_success) return
            if (j==1) then
                nf = size(one,1)
                nt = size(one,2)
                deallocate(values)
                allocate(values(nf,state%nchannels,nt))
                times = t
            else
                if (size(one,1) /= nf .or. size(one,2) /= nt) then
                    status = spectral_status_invalid_input
                    return
                end if
            end if
            values(:,j,:) = one
        end do
        state%finished = .true.
    end subroutine forward_finish_core

    subroutine stft_stream_multi_finish(state,values,status,times)
        type(stft_stream_multi), intent(inout) :: state !! active forward stream
        complex(dp), allocatable, intent(out) :: values(:,:,:) !! final spectral frames
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! shared frame coordinates
        real(dp), allocatable :: t(:)

        call forward_finish_core(state,values,t,status,.false.)
        if (present(times)) times = t
    end subroutine stft_stream_multi_finish

    subroutine spectrogram_stream_multi_finish(state,power,status,times)
        type(stft_stream_multi), intent(inout) :: state !! active spectrogram stream
        real(dp), allocatable, intent(out) :: power(:,:,:) !! final PSD frames
        integer, intent(out) :: status !! spectral status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! shared frame coordinates
        complex(dp), allocatable :: values(:,:,:)
        real(dp), allocatable :: t(:)

        call forward_finish_core(state,values,t,status,.true.)
        power = real(values,dp)
        if (present(times)) times = t
    end subroutine spectrogram_stream_multi_finish

    subroutine istft_stream_multi_init(state,nchannels,nfreq,status,fs,window,nperseg, &
        noverlap,nfft,boundary,scaling,window_values,input_onesided)
        type(istft_stream_multi), intent(inout) :: state !! unchanged on invalid configuration
        integer, intent(in) :: nchannels !! number of independent output signals
        integer, intent(in) :: nfreq !! number of input STFT bins
        integer, intent(out) :: status !! inverse STFT status code
        real(dp), optional, intent(in) :: fs !! sampling frequency
        character(*), optional, intent(in) :: window !! synthesis window
        integer, optional, intent(in) :: nperseg !! synthesis window length
        integer, optional, intent(in) :: noverlap !! overlap samples
        integer, optional, intent(in) :: nfft !! transform length
        logical, optional, intent(in) :: boundary !! trim analysis boundary padding
        character(*), optional, intent(in) :: scaling !! spectrum or psd
        real(dp), optional, intent(in) :: window_values(:) !! synthesis window values
        logical, optional, intent(in) :: input_onesided !! real-signal positive-frequency bins
        type(istft_stream_multi) :: trial
        integer :: j

        status = istft_status_invalid_input
        if (nchannels < 1) return
        allocate(trial%channels(nchannels))
        do j=1,nchannels
            call istft_stream_init(trial%channels(j),nfreq,status,fs,window,nperseg, &
                noverlap,nfft,boundary,scaling,window_values,input_onesided)
            if (status /= istft_status_success) return
        end do
        trial%nchannels = nchannels
        trial%nfreq = nfreq
        trial%onesided = .true.
        if (present(input_onesided)) trial%onesided = input_onesided
        trial%active = .true.
        state = trial
    end subroutine istft_stream_multi_init

    subroutine inverse_push_core(state,frame,samples,times,status,real_output)
        type(istft_stream_multi), intent(inout) :: state !! active inverse stream
        complex(dp), intent(in) :: frame(:,:) !! frequency by channel input frame
        complex(dp), allocatable, intent(out) :: samples(:,:) !! completed samples by channel
        real(dp), allocatable, intent(out) :: times(:) !! common sample coordinates
        integer, intent(out) :: status !! inverse STFT status code
        logical, intent(in) :: real_output !! request real-valued reconstruction
        complex(dp), allocatable :: one(:)
        real(dp), allocatable :: r(:), t(:)
        integer :: j, ns, st
        logical :: noninvertible

        allocate(samples(0,0),times(0))
        status = istft_status_invalid_input
        if (.not.state%active .or. state%finished) return
        if (real_output .and. .not.state%onesided) return
        if (size(frame,1) /= state%nfreq .or. size(frame,2) /= state%nchannels) return
        if (.not.all(ieee_is_finite(real(frame,dp)))) return
        if (.not.all(ieee_is_finite(aimag(frame)))) return
        noninvertible = .false.
        do j=1,state%nchannels
            if (real_output) then
                call istft_stream_push(state%channels(j),frame(:,j),r,st,t)
                one = cmplx(r,0.0_dp,dp)
            else
                call istft_stream_push(state%channels(j),frame(:,j),one,st,t)
            end if
            if (st == istft_status_invalid_input) return
            noninvertible = noninvertible .or. st==istft_status_noninvertible
            if (j==1) then
                ns = size(one)
                deallocate(samples)
                allocate(samples(ns,state%nchannels))
                times = t
            else
                if (size(one) /= ns) return
            end if
            samples(:,j) = one
        end do
        status = istft_status_success
        if (noninvertible) status = istft_status_noninvertible
    end subroutine inverse_push_core

    subroutine inverse_push_real(state,frame,samples,status,times)
        type(istft_stream_multi), intent(inout) :: state !! inverse stream
        complex(dp), intent(in) :: frame(:,:) !! frequency by channel frame
        real(dp), allocatable, intent(out) :: samples(:,:) !! completed real samples
        integer, intent(out) :: status !! inverse STFT status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! sample coordinates
        complex(dp), allocatable :: values(:,:)
        real(dp), allocatable :: t(:)

        call inverse_push_core(state,frame,values,t,status,.true.)
        samples = real(values,dp)
        if (present(times)) times = t
    end subroutine inverse_push_real

    subroutine inverse_push_complex(state,frame,samples,status,times)
        type(istft_stream_multi), intent(inout) :: state !! inverse stream
        complex(dp), intent(in) :: frame(:,:) !! frequency by channel frame
        complex(dp), allocatable, intent(out) :: samples(:,:) !! completed complex samples
        integer, intent(out) :: status !! inverse STFT status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! sample coordinates
        real(dp), allocatable :: t(:)

        call inverse_push_core(state,frame,samples,t,status,.false.)
        if (present(times)) times = t
    end subroutine inverse_push_complex

    subroutine inverse_finish_core(state,samples,times,status,real_output)
        type(istft_stream_multi), intent(inout) :: state !! inverse stream
        complex(dp), allocatable, intent(out) :: samples(:,:) !! remaining completed samples
        real(dp), allocatable, intent(out) :: times(:) !! sample coordinates
        integer, intent(out) :: status !! inverse STFT status code
        logical, intent(in) :: real_output !! request real-valued samples
        complex(dp), allocatable :: one(:)
        real(dp), allocatable :: r(:), t(:)
        integer :: j, ns, st
        logical :: noninvertible

        allocate(samples(0,0),times(0))
        status = istft_status_invalid_input
        if (.not.state%active .or. state%finished) return
        if (real_output .and. .not.state%onesided) return
        noninvertible = .false.
        do j=1,state%nchannels
            if (real_output) then
                call istft_stream_finish(state%channels(j),r,st,t)
                one = cmplx(r,0.0_dp,dp)
            else
                call istft_stream_finish(state%channels(j),one,st,t)
            end if
            if (st == istft_status_invalid_input) return
            noninvertible = noninvertible .or. st==istft_status_noninvertible
            if (j==1) then
                ns = size(one)
                deallocate(samples)
                allocate(samples(ns,state%nchannels))
                times = t
            else
                if (size(one) /= ns) return
            end if
            samples(:,j) = one
        end do
        state%finished = .true.
        status = istft_status_success
        if (noninvertible) status = istft_status_noninvertible
    end subroutine inverse_finish_core

    subroutine inverse_finish_real(state,samples,status,times)
        type(istft_stream_multi), intent(inout) :: state !! inverse stream
        real(dp), allocatable, intent(out) :: samples(:,:) !! completed samples by channel
        integer, intent(out) :: status !! inverse STFT status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! sample coordinates
        complex(dp), allocatable :: values(:,:)
        real(dp), allocatable :: t(:)

        call inverse_finish_core(state,values,t,status,.true.)
        samples = real(values,dp)
        if (present(times)) times = t
    end subroutine inverse_finish_real

    subroutine inverse_finish_complex(state,samples,status,times)
        type(istft_stream_multi), intent(inout) :: state !! inverse stream
        complex(dp), allocatable, intent(out) :: samples(:,:) !! completed samples by channel
        integer, intent(out) :: status !! inverse STFT status code
        real(dp), allocatable, optional, intent(out) :: times(:) !! sample coordinates
        real(dp), allocatable :: t(:)

        call inverse_finish_core(state,samples,t,status,.false.)
        if (present(times)) times = t
    end subroutine inverse_finish_complex

end module scifort_signal_stream_multi
