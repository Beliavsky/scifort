! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style inverse short-time Fourier transform and overlap-add checks.
! Independent implementation of weighted overlap-add reconstruction.
module scifort_signal_istft
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_named_window, forward_fft, array_median, &
        spectral_status_success
    implicit none
    private

    integer, parameter, public :: istft_status_success = 0
    integer, parameter, public :: istft_status_invalid_input = 1
    integer, parameter, public :: istft_status_noninvertible = 2
    real(dp), parameter :: norm_tolerance = 1.0e-10_dp

    interface istft
        module procedure istft_real
        module procedure istft_complex
    end interface
    interface check_cola
        module procedure check_cola_named
        module procedure check_cola_values
    end interface
    interface check_nola
        module procedure check_nola_named
        module procedure check_nola_values
    end interface
    public :: istft, check_cola, check_nola

contains

    subroutine istft_real(zxx, samples, status, fs, window, nperseg, noverlap, nfft, &
            boundary, scaling, times, window_values, input_onesided, &
            window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: zxx(:,:) !! STFT, frequency by frame
        real(dp), allocatable, intent(out) :: samples(:) !! reconstructed real signal
        integer, intent(out) :: status !! success, invalid input, or incomplete overlap
        real(dp), intent(in), optional :: fs !! samples per unit time; default 1
        character(*), intent(in), optional :: window !! named periodic window; default hann
        integer, intent(in), optional :: nperseg !! window size; inferred when omitted
        integer, intent(in), optional :: noverlap !! overlapped samples; default half window
        integer, intent(in), optional :: nfft !! FFT size; inferred when omitted
        logical, intent(in), optional :: boundary !! remove window-center extensions; default true
        character(*), intent(in), optional :: scaling !! spectrum (default) or psd
        real(dp), allocatable, intent(out), optional :: times(:) !! zero-based sample times
        real(dp), intent(in), optional :: window_values(:) !! custom synthesis/analysis window
        logical, intent(in), optional :: input_onesided !! must be true for real output
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        complex(dp), allocatable :: reconstructed(:)
        real(dp), allocatable :: time_grid(:)
        logical :: onesided

        allocate(samples(0))
        if (present(times)) allocate(times(0))
        onesided = .true.
        if (present(input_onesided)) onesided = input_onesided
        status = istft_status_invalid_input
        if (.not. onesided) return
        call inverse_stft_kernel(zxx, reconstructed, time_grid, status, &
            fs, window, nperseg, noverlap, nfft, onesided, boundary, scaling, window_values, &
            window_parameter, window_parameter2, window_coefficients)
        if (status == istft_status_invalid_input) return
        samples = real(reconstructed, dp)
        if (present(times)) times = time_grid
    end subroutine istft_real

    subroutine istft_complex(zxx, samples, status, fs, window, nperseg, noverlap, nfft, &
            boundary, scaling, times, window_values, input_onesided, &
            window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: zxx(:,:) !! STFT, frequency by frame
        complex(dp), allocatable, intent(out) :: samples(:) !! reconstructed complex signal
        integer, intent(out) :: status !! success, invalid input, or incomplete overlap
        real(dp), intent(in), optional :: fs !! samples per unit time; default 1
        character(*), intent(in), optional :: window !! named periodic window; default hann
        integer, intent(in), optional :: nperseg !! window size; inferred when omitted
        integer, intent(in), optional :: noverlap !! overlapped samples; default half window
        integer, intent(in), optional :: nfft !! FFT size; inferred when omitted
        logical, intent(in), optional :: boundary !! remove window-center extensions; default true
        character(*), intent(in), optional :: scaling !! spectrum (default) or psd
        real(dp), allocatable, intent(out), optional :: times(:) !! zero-based sample times
        real(dp), intent(in), optional :: window_values(:) !! custom synthesis/analysis window
        logical, intent(in), optional :: input_onesided !! default false for complex output
        real(dp), intent(in), optional :: window_parameter !! named window shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second named window parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), allocatable :: time_grid(:)
        logical :: onesided

        onesided = .false.
        if (present(input_onesided)) onesided = input_onesided
        call inverse_stft_kernel(zxx, samples, time_grid, status, &
            fs, window, nperseg, noverlap, nfft, onesided, boundary, scaling, window_values, &
            window_parameter, window_parameter2, window_coefficients)
        if (present(times)) times = time_grid
    end subroutine istft_complex

    subroutine inverse_stft_kernel(zxx, samples, times, status, fs, window, nperseg, &
            noverlap, nfft, onesided, boundary, scaling, window_values, &
            window_parameter, window_parameter2, window_coefficients)
        complex(dp), intent(in) :: zxx(:,:) !! frequency-bin by frame array
        complex(dp), allocatable, intent(out) :: samples(:) !! reconstructed complex samples
        real(dp), allocatable, intent(out) :: times(:) !! corresponding zero-based time points
        integer, intent(out) :: status !! inverse-transform status
        real(dp), intent(in), optional :: fs !! sampling frequency
        character(*), intent(in), optional :: window !! periodic window name
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap length
        integer, intent(in), optional :: nfft !! full transform length
        logical, intent(in) :: onesided !! use Hermitian symmetry
        logical, intent(in), optional :: boundary !! discard extensions at both ends
        character(*), intent(in), optional :: scaling !! spectrum or psd
        real(dp), intent(in), optional :: window_values(:) !! explicit window
        real(dp), intent(in), optional :: window_parameter !! window family shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! second shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! cosine series weights

        integer :: nfreq, nframes, default_nfft, nw, nf, no, hop, nall, edge_pad, nout
        integer :: i, j, k, start, ierr
        real(dp) :: rate, factor
        real(dp), allocatable :: w(:), normalizer(:)
        complex(dp), allocatable :: bins(:), output(:)
        character(:), allocatable :: win_name, scale
        logical :: edges

        allocate(samples(0), times(0))
        status = istft_status_invalid_input
        nfreq = size(zxx, 1)
        nframes = size(zxx, 2)
        if (nfreq < 1 .or. nframes < 1) return
        if (.not. all(ieee_is_finite(real(zxx,dp)))) return
        if (.not. all(ieee_is_finite(aimag(zxx)))) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate) .or. rate <= 0.0_dp) return

        default_nfft = nfreq
        if (onesided) default_nfft = 2*(nfreq - 1)
        if (onesided .and. nfreq == 1) default_nfft = 1
        nw = default_nfft
        if (present(nperseg)) nw = nperseg
        if (nw < 1) return
        nf = default_nfft
        if (onesided .and. nw == default_nfft + 1) nf = nw
        if (present(nfft)) nf = nfft
        if (nf < nw .or. nf < 1) return
        if (onesided) then
            if (nfreq /= nf / 2 + 1) return
        else
            if (nfreq /= nf) return
        end if
        no = nw / 2
        if (present(noverlap)) no = noverlap
        if (no < 0 .or. no >= nw) return
        hop = nw - no
        edges = .true.
        if (present(boundary)) edges = boundary
        scale = 'spectrum'
        if (present(scaling)) scale = trim(scaling)
        if (scale /= 'spectrum' .and. scale /= 'psd') return
        if (present(window_values)) then
            if (present(window) .or. present(window_parameter) .or. &
                present(window_parameter2) .or. present(window_coefficients)) return
        end if

        if (present(window_values)) then
            if (size(window_values) /= nw) return
            w = window_values
        else
            win_name = 'hann'
            if (present(window)) win_name = trim(window)
            w = spectral_named_window(win_name, nw, ierr, window_parameter, &
                window_parameter2, window_coefficients)
            if (ierr /= spectral_status_success) return
        end if
        if (.not. all(ieee_is_finite(w))) return
        if (scale == 'spectrum') then
            factor = sum(w)
        else
            factor = sqrt(rate*sum(w*w))
        end if
        if (.not. ieee_is_finite(factor)) return
        ! On normal SciPy round trips these are positive. A zero window is
        ! technically processable, but its overlap-add normalization vanishes.
        if (nframes-1 > (huge(nframes)-nw) / hop) return
        nall = nw + (nframes - 1)*hop
        edge_pad = 0
        if (edges) edge_pad = nw/2
        nout = nall - 2*edge_pad
        if (nout < 0) return
        allocate(output(nall), normalizer(nall), bins(nf))
        output = (0.0_dp,0.0_dp)
        normalizer = 0.0_dp
        do j = 1, nframes
            bins = (0.0_dp,0.0_dp)
            if (onesided) then
                bins(1:nfreq) = zxx(:,j)
                do k = 2, nfreq
                    i = nf - k + 2
                    if (i /= k) bins(i) = conjg(zxx(k,j))
                end do
            else
                bins = zxx(:,j)
            end if
            ! conjg(FFT(conjg(x))) / nf is the inverse for all FFT lengths.
            bins = conjg(bins)
            call forward_fft(bins)
            bins = conjg(bins) / real(nf,dp)
            start = (j-1)*hop
            do i = 1, nw
                output(start+i) = output(start+i) + bins(i)*factor*w(i)
                normalizer(start+i) = normalizer(start+i) + w(i)*w(i)
            end do
        end do
        deallocate(samples, times)
        allocate(samples(nout), times(nout))
        status = istft_status_success
        do i = 1, nout
            k = i + edge_pad
            if (normalizer(k) > norm_tolerance) then
                samples(i) = output(k) / normalizer(k)
            else
                samples(i) = output(k)
                status = istft_status_noninvertible
            end if
            times(i) = real(i-1, dp) / rate
        end do
        if (onesided) samples = cmplx(real(samples,dp),0.0_dp,dp)
    end subroutine inverse_stft_kernel

    logical function check_cola_named(window, nperseg, noverlap, tol, status, &
            window_parameter, window_parameter2, window_coefficients) result(valid)
        character(*), intent(in) :: window !! named periodic window
        integer, intent(in) :: nperseg !! segment length
        integer, intent(in) :: noverlap !! overlap length
        real(dp), intent(in), optional :: tol !! absolute tolerance, default 1e-10
        integer, intent(out), optional :: status !! zero on valid input
        real(dp), intent(in), optional :: window_parameter !! named window first parameter
        real(dp), intent(in), optional :: window_parameter2 !! named window second parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! cosine family coefficients
        real(dp), allocatable :: w(:)
        integer :: ierr

        valid = .false.
        if (present(status)) status = istft_status_invalid_input
        if (nperseg < 1 .or. noverlap < 0 .or. noverlap >= nperseg) return
        w = spectral_named_window(window, nperseg, ierr, window_parameter, &
            window_parameter2, window_coefficients)
        if (ierr /= spectral_status_success) return
        valid = overlap_condition(w, noverlap, .false., tol)
        if (present(status)) status = istft_status_success
    end function check_cola_named

    logical function check_nola_named(window, nperseg, noverlap, tol, status, &
            window_parameter, window_parameter2, window_coefficients) result(valid)
        character(*), intent(in) :: window !! named periodic window
        integer, intent(in) :: nperseg !! segment length
        integer, intent(in) :: noverlap !! overlap length
        real(dp), intent(in), optional :: tol !! absolute tolerance, default 1e-10
        integer, intent(out), optional :: status !! zero on valid input
        real(dp), intent(in), optional :: window_parameter !! named window first parameter
        real(dp), intent(in), optional :: window_parameter2 !! named window second parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! cosine family coefficients
        real(dp), allocatable :: w(:)
        integer :: ierr

        valid = .false.
        if (present(status)) status = istft_status_invalid_input
        if (nperseg < 1 .or. noverlap < 0 .or. noverlap >= nperseg) return
        w = spectral_named_window(window, nperseg, ierr, window_parameter, &
            window_parameter2, window_coefficients)
        if (ierr /= spectral_status_success) return
        valid = overlap_condition(w, noverlap, .true., tol)
        if (present(status)) status = istft_status_success
    end function check_nola_named

    logical function check_cola_values(window, nperseg, noverlap, tol, status) result(valid)
        real(dp), intent(in) :: window(:) !! explicit periodic or symmetric window samples
        integer, intent(in) :: nperseg !! segment length
        integer, intent(in) :: noverlap !! overlap length
        real(dp), intent(in), optional :: tol !! absolute tolerance, default 1e-10
        integer, intent(out), optional :: status !! zero on valid input

        valid = .false.
        if (present(status)) status = istft_status_invalid_input
        if (size(window) /= nperseg .or. nperseg < 1) return
        if (noverlap < 0 .or. noverlap >= nperseg) return
        if (.not. all(ieee_is_finite(window))) return
        valid = overlap_condition(window, noverlap, .false., tol)
        if (present(status)) status = istft_status_success
    end function check_cola_values

    logical function check_nola_values(window, nperseg, noverlap, tol, status) result(valid)
        real(dp), intent(in) :: window(:) !! explicit periodic or symmetric window samples
        integer, intent(in) :: nperseg !! segment length
        integer, intent(in) :: noverlap !! overlap length
        real(dp), intent(in), optional :: tol !! absolute tolerance, default 1e-10
        integer, intent(out), optional :: status !! zero on valid input

        valid = .false.
        if (present(status)) status = istft_status_invalid_input
        if (size(window) /= nperseg .or. nperseg < 1) return
        if (noverlap < 0 .or. noverlap >= nperseg) return
        if (.not. all(ieee_is_finite(window))) return
        valid = overlap_condition(window, noverlap, .true., tol)
        if (present(status)) status = istft_status_success
    end function check_nola_values

    logical function overlap_condition(w, noverlap, squared, tol) result(valid)
        real(dp), intent(in) :: w(:) !! window samples
        integer, intent(in) :: noverlap !! overlap length
        logical, intent(in) :: squared !! NOLA uses squares; COLA uses values
        real(dp), intent(in), optional :: tol !! absolute numerical tolerance
        real(dp), allocatable :: sums(:)
        real(dp) :: threshold, med
        integer :: i, hop

        threshold = norm_tolerance
        if (present(tol)) threshold = tol
        valid = .false.
        if (.not. ieee_is_finite(threshold)) return
        hop = size(w) - noverlap
        allocate(sums(hop))
        sums = 0.0_dp
        do i = 1, size(w)
            if (squared) then
                sums(mod(i-1,hop)+1) = sums(mod(i-1,hop)+1) + w(i)*w(i)
            else
                sums(mod(i-1,hop)+1) = sums(mod(i-1,hop)+1) + w(i)
            end if
        end do
        if (squared) then
            valid = minval(sums) > threshold
        else
            med = array_median(sums)
            valid = maxval(abs(sums-med)) < threshold
        end if
    end function overlap_condition

end module scifort_signal_istft
