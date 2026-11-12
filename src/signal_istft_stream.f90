! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Incremental inverse STFT for scipy.signal.stft-style frames.
! Retains one window and the pending overlap, not the full reconstruction.
module scifort_signal_istft_stream
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan, get_window, forward_fft_planned, spectral_status_success
    use scifort_signal_istft, only : istft_status_success, istft_status_invalid_input, &
        istft_status_noninvertible
    implicit none
    private

    real(dp), parameter :: norm_tolerance = 1.0e-10_dp

    type, public :: istft_stream
        private
        integer :: nfreq = 0, nperseg = 0, nfft = 0, hop = 0
        integer :: edge_pad = 0, ring_size = 0
        integer :: frames = 0, consumed = 0, produced = 0
        real(dp) :: fs = 1.0_dp, factor = 1.0_dp
        logical :: onesided = .true., active = .false., finished = .false.
        logical :: noninvertible = .false.
        real(dp), allocatable :: w(:), normalizer(:)
        complex(dp), allocatable :: overlap(:), bins(:)
        type(fft_plan) :: transform_plan
    end type istft_stream

    interface istft_stream_push
        module procedure stream_push_real, stream_push_complex
    end interface
    interface istft_stream_finish
        module procedure stream_finish_real, stream_finish_complex
    end interface
    public :: istft_stream_init, istft_stream_push, istft_stream_finish

contains

    subroutine istft_stream_init(state, nfreq, status, fs, window, nperseg, noverlap, &
        nfft, boundary, scaling, window_values, input_onesided)
        type(istft_stream), intent(inout) :: state !! streaming reconstruction state
        integer, intent(in) :: nfreq !! number of input STFT frequency bins
        integer, intent(out) :: status !! invalid parameters or success
        real(dp), intent(in), optional :: fs !! sampling frequency, default 1
        character(*), intent(in), optional :: window !! periodic window name, default hann
        integer, intent(in), optional :: nperseg !! length of the synthesis window
        integer, intent(in), optional :: noverlap !! overlapping samples
        integer, intent(in), optional :: nfft !! FFT length, inferred if absent
        logical, intent(in), optional :: boundary !! discard padded edges, default true
        character(*), intent(in), optional :: scaling !! spectrum or psd
        real(dp), intent(in), optional :: window_values(:) !! explicit synthesis window
        logical, intent(in), optional :: input_onesided !! one-sided input, default true
        type(istft_stream) :: trial
        integer :: default_nfft, nw, nf, no, ierr
        real(dp) :: rate, factor
        logical :: onesided, edges
        character(:), allocatable :: win_name, scale
        real(dp), allocatable :: w(:)

        status = istft_status_invalid_input
        if (nfreq < 1) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate) .or. rate <= 0.0_dp) return
        onesided = .true.
        if (present(input_onesided)) onesided = input_onesided
        default_nfft = nfreq
        if (onesided) then
            if (nfreq > huge(nfreq)/2) return
            default_nfft = max(1,2*(nfreq-1))
        end if
        nw = default_nfft
        if (present(nperseg)) nw = nperseg
        if (nw < 1) return
        nf = default_nfft
        if (onesided .and. default_nfft < huge(nw)) then
            if (nw == default_nfft+1) nf = nw
        end if
        if (present(nfft)) nf = nfft
        if (nf < nw .or. nf < 1) return
        if (onesided) then
            if (nfreq /= nf/2+1) return
        else
            if (nfreq /= nf) return
        end if
        no = nw/2
        if (present(noverlap)) no = noverlap
        if (no < 0 .or. no >= nw) return
        edges = .true.
        if (present(boundary)) edges = boundary
        scale = 'spectrum'
        if (present(scaling)) scale = trim(scaling)
        if (scale /= 'spectrum' .and. scale /= 'psd') return
        if (present(window) .and. present(window_values)) return
        if (present(window_values)) then
            if (size(window_values) /= nw) return
            w = window_values
        else
            win_name = 'hann'
            if (present(window)) win_name = trim(window)
            w = get_window(win_name, nw, status=ierr)
            if (ierr /= spectral_status_success) return
        end if
        if (.not. all(ieee_is_finite(w))) return
        if (scale == 'spectrum') then
            factor = sum(w)
        else
            factor = sqrt(rate*sum(w*w))
        end if
        if (.not. ieee_is_finite(factor)) return

        trial%nfreq = nfreq
        trial%nfft = nf
        trial%nperseg = nw
        trial%hop = nw - no
        trial%edge_pad = 0
        if (edges) trial%edge_pad = nw/2
        if (nw > huge(nw)-trial%edge_pad) return
        trial%ring_size = nw + trial%edge_pad
        trial%fs = rate
        trial%factor = factor
        trial%onesided = onesided
        trial%w = w
        allocate(trial%overlap(trial%ring_size), trial%normalizer(trial%ring_size))
        allocate(trial%bins(nf))
        trial%overlap = cmplx(0.0_dp,0.0_dp,dp)
        trial%normalizer = 0.0_dp
        trial%active = .true.
        call trial%transform_plan%initialize(nf, ierr)
        if (ierr /= spectral_status_success) return
        state = trial
        status = istft_status_success
    end subroutine istft_stream_init

    subroutine stream_frame(state, frame, samples, times, status)
        type(istft_stream), intent(inout) :: state !! previously initialized state
        complex(dp), intent(in) :: frame(:) !! one frequency-domain frame
        complex(dp), allocatable, intent(out) :: samples(:) !! completed reconstructed samples
        real(dp), allocatable, intent(out) :: times(:) !! completed sample time coordinates
        integer, intent(out) :: status !! success, noninvertible, or invalid input
        integer :: k, idx, start, emit_end, j, i, n

        allocate(samples(0), times(0))
        status = istft_status_invalid_input
        if (.not. state%active .or. state%finished) return
        if (size(frame) /= state%nfreq) return
        if (.not. all(ieee_is_finite(real(frame,dp)))) return
        if (.not. all(ieee_is_finite(aimag(frame)))) return
        if (state%frames > (huge(start)-state%nperseg)/state%hop) return

        state%bins = cmplx(0.0_dp,0.0_dp,dp)
        if (state%onesided) then
            state%bins(1:state%nfreq) = frame
            do k = 2,state%nfreq
                j = state%nfft-k+2
                if (j /= k) state%bins(j) = conjg(frame(k))
            end do
        else
            state%bins = frame
        end if
        state%bins = conjg(state%bins)
        call forward_fft_planned(state%bins, state%transform_plan)
        state%bins = conjg(state%bins) / real(state%nfft,dp)
        start = state%frames*state%hop
        do i = 1,state%nperseg
            idx = mod(start+i-1,state%ring_size)+1
            state%overlap(idx) = state%overlap(idx) + state%bins(i)*state%factor*state%w(i)
            state%normalizer(idx) = state%normalizer(idx) + state%w(i)*state%w(i)
        end do
        state%frames = state%frames+1
        emit_end = max(0,state%frames*state%hop-state%edge_pad)
        call collect(state,emit_end,emit_end,samples,times)
        status = istft_status_success
        if (state%noninvertible) status = istft_status_noninvertible
    end subroutine stream_frame

    subroutine collect(state, clear_end, keep_end, samples, times)
        type(istft_stream), intent(inout) :: state !! current overlap and normalization state
        integer, intent(in) :: clear_end !! raw samples which are now final, exclusive
        integer, intent(in) :: keep_end !! raw samples before this bound are retained
        complex(dp), allocatable, intent(out) :: samples(:) !! completed samples
        real(dp), allocatable, intent(out) :: times(:) !! completed times
        integer :: first, n, raw, index, j, k
        complex(dp) :: z

        first = max(state%consumed,state%edge_pad)
        n = max(0,keep_end-first)
        allocate(samples(n),times(n))
        j = 0
        do raw = state%consumed,clear_end-1
            index = mod(raw,state%ring_size)+1
            if (raw >= state%edge_pad .and. raw < keep_end) then
                z = state%overlap(index)
                if (state%normalizer(index) > norm_tolerance) then
                    z = z / state%normalizer(index)
                else
                    state%noninvertible = .true.
                end if
                if (state%onesided) z = cmplx(real(z,dp),0.0_dp,dp)
                j = j+1
                samples(j) = z
                k = state%produced+j-1
                times(j) = real(k,dp)/state%fs
            end if
            state%overlap(index) = cmplx(0.0_dp,0.0_dp,dp)
            state%normalizer(index) = 0.0_dp
        end do
        state%consumed = clear_end
        state%produced = state%produced+n
    end subroutine collect

    subroutine stream_push_complex(state, frame, samples, status, times)
        type(istft_stream), intent(inout) :: state !! inverse STFT stream
        complex(dp), intent(in) :: frame(:) !! complex STFT frame
        complex(dp), allocatable, intent(out) :: samples(:) !! finalized samples
        integer, intent(out) :: status !! status of processing so far
        real(dp), allocatable, intent(out), optional :: times(:) !! corresponding time values
        real(dp), allocatable :: chunk_times(:)

        call stream_frame(state,frame,samples,chunk_times,status)
        if (present(times)) times = chunk_times
    end subroutine stream_push_complex

    subroutine stream_push_real(state, frame, samples, status, times)
        type(istft_stream), intent(inout) :: state !! inverse STFT stream
        complex(dp), intent(in) :: frame(:) !! one-sided STFT frame
        real(dp), allocatable, intent(out) :: samples(:) !! finalized real samples
        integer, intent(out) :: status !! status of processing so far
        real(dp), allocatable, intent(out), optional :: times(:) !! corresponding time values
        complex(dp), allocatable :: chunk(:)
        real(dp), allocatable :: chunk_times(:)

        allocate(samples(0))
        if (.not. state%onesided) then
            status = istft_status_invalid_input
            if (present(times)) allocate(times(0))
            return
        end if
        call stream_frame(state,frame,chunk,chunk_times,status)
        samples = real(chunk,dp)
        if (present(times)) times = chunk_times
    end subroutine stream_push_real

    subroutine stream_end(state, samples, times, status)
        type(istft_stream), intent(inout) :: state !! inverse STFT stream
        complex(dp), allocatable, intent(out) :: samples(:) !! remaining reconstructed samples
        real(dp), allocatable, intent(out) :: times(:) !! corresponding time values
        integer, intent(out) :: status !! success, noninvertible, or invalid input
        integer :: nall, trim_end

        allocate(samples(0),times(0))
        status = istft_status_invalid_input
        if (.not. state%active .or. state%finished .or. state%frames < 1) return
        nall = state%nperseg+(state%frames-1)*state%hop
        trim_end = nall-state%edge_pad
        call collect(state,nall,trim_end,samples,times)
        state%finished = .true.
        status = istft_status_success
        if (state%noninvertible) status = istft_status_noninvertible
    end subroutine stream_end

    subroutine stream_finish_complex(state, samples, status, times)
        type(istft_stream), intent(inout) :: state !! inverse STFT stream
        complex(dp), allocatable, intent(out) :: samples(:) !! remaining complex samples
        integer, intent(out) :: status !! final stream status
        real(dp), allocatable, intent(out), optional :: times(:) !! corresponding time values
        real(dp), allocatable :: chunk_times(:)

        call stream_end(state,samples,chunk_times,status)
        if (present(times)) times = chunk_times
    end subroutine stream_finish_complex

    subroutine stream_finish_real(state, samples, status, times)
        type(istft_stream), intent(inout) :: state !! inverse STFT stream
        real(dp), allocatable, intent(out) :: samples(:) !! remaining real samples
        integer, intent(out) :: status !! final stream status
        real(dp), allocatable, intent(out), optional :: times(:) !! corresponding time values
        complex(dp), allocatable :: chunk(:)
        real(dp), allocatable :: chunk_times(:)

        allocate(samples(0))
        if (.not. state%onesided) then
            status = istft_status_invalid_input
            if (present(times)) allocate(times(0))
            return
        end if
        call stream_end(state,chunk,chunk_times,status)
        samples = real(chunk,dp)
        if (present(times)) times = chunk_times
    end subroutine stream_finish_real

end module scifort_signal_istft_stream
