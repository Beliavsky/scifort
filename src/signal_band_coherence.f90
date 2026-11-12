! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Exact piecewise-linear integration of spectral matrices over frequency bands.
! Band coherence is computed from the integrated CSD, not the average of
! frequency-by-frequency coherence; this is a SciFort extension, not SciPy API.
module scifort_signal_band_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_result, &
        partial_coherence_tf_result, partial_coherence_matrix, partial_coherence_tf_matrix
    use scifort_signal_multiple_coherence, only : multiple_coherence_result, &
        multiple_coherence_tf_result, multiple_coherence_matrix, multiple_coherence_tf_matrix
    implicit none
    private

    interface multitaper_band_multiple_coherence
        module procedure multitaper_band_multiple_coherence_real
        module procedure multitaper_band_multiple_coherence_complex
    end interface
    interface multitaper_band_partial_coherence
        module procedure multitaper_band_partial_coherence_real
        module procedure multitaper_band_partial_coherence_complex
    end interface
    interface multitaper_tf_band_partial_coherence
        module procedure multitaper_tf_band_partial_coherence_real
        module procedure multitaper_tf_band_partial_coherence_complex
    end interface
    interface multitaper_tf_band_multiple_coherence
        module procedure multitaper_tf_band_multiple_coherence_real
        module procedure multitaper_tf_band_multiple_coherence_complex
    end interface

    public :: integrate_spectral_bands, integrate_tf_spectral_bands
    public :: band_multiple_coherence_matrix, band_multiple_coherence_tf_matrix
    public :: band_partial_coherence_matrix, band_partial_coherence_tf_matrix
    public :: multitaper_band_multiple_coherence, multitaper_tf_band_multiple_coherence
    public :: multitaper_band_partial_coherence, multitaper_tf_band_partial_coherence

contains

    ! Stable O(n log n) merge sort. This also accepts the FFT ordering used
    ! for complex signals: [zero, positive frequencies, negative frequencies].
    pure subroutine frequency_order(frequencies, order, ok)
        real(dp), intent(in) :: frequencies(:) !! finite input frequency grid, in Hz
        integer, intent(out) :: order(:) !! indices that sort the grid ascending
        logical, intent(out) :: ok !! true for finite, distinct frequencies
        integer :: work(size(order)), n, width, left, middle, right, i, j, k

        n = size(frequencies)
        ok = .false.
        if (size(order) /= n .or. n < 2) return
        if (.not. all(ieee_is_finite(frequencies))) return
        order = [(i,i=1,n)]
        width = 1
        do while (width < n)
            left = 1
            do while (left <= n)
                middle = min(n,left+width-1)
                right = min(n,middle+width)
                i = left
                j = middle+1
                do k = left, right
                    if (i > middle) then
                        work(k) = order(j)
                        j = j+1
                    else if (j > right) then
                        work(k) = order(i)
                        i = i+1
                    else if (frequencies(order(i)) <= frequencies(order(j))) then
                        work(k) = order(i)
                        i = i+1
                    else
                        work(k) = order(j)
                        j = j+1
                    end if
                end do
                left = right+1
            end do
            order = work
            if (width > n/2) exit
            width = 2*width
        end do
        do k = 2, n
            if (frequencies(order(k)) <= frequencies(order(k-1))) return
        end do
        ok = .true.
    end subroutine frequency_order

    ! Integrate the linearly interpolated CSD over each requested band.
    ! For each clipped interval [a,b] of [f_i,f_{i+1}], its two trapezoidal
    ! interpolation weights sum to b-a. No interpolated spectrum is evaluated
    ! outside the supplied frequency grid; endpoints are included exactly.
    pure subroutine integration_weights(frequencies, bands, weights, ok)
        real(dp), intent(in) :: frequencies(:) !! FFT frequency coordinates, any permutation
        real(dp), intent(in) :: bands(:,:) !! requested limits (2,nband), in Hz
        real(dp), intent(out) :: weights(:,:) !! quadrature weights (frequency,band)
        logical, intent(out) :: ok !! valid grids and in-range nonempty bands
        integer :: order(size(frequencies)), nf, nb, k, ib, i, j
        real(dp) :: f0, f1, a, b, h, mean_fraction, length

        ok = .false.
        weights = 0.0_dp
        nf = size(frequencies)
        if (size(bands,1) /= 2) return
        nb = size(bands,2)
        if (nb < 1) return
        if (size(weights,1) /= nf .or. size(weights,2) /= nb) return
        if (.not. all(ieee_is_finite(bands))) return
        call frequency_order(frequencies,order,ok)
        if (.not. ok) return
        ok = .false.
        do ib = 1, nb
            if (bands(1,ib) >= bands(2,ib)) return
            if (bands(1,ib) < frequencies(order(1))) return
            if (bands(2,ib) > frequencies(order(nf))) return
        end do
        do ib = 1, nb
            do k = 1, nf-1
                i = order(k)
                j = order(k+1)
                f0 = frequencies(i)
                f1 = frequencies(j)
                a = max(f0,bands(1,ib))
                b = min(f1,bands(2,ib))
                if (b <= a) cycle
                h = f1-f0
                length = b-a
                ! Average local normalized abscissa of a linear function.
                mean_fraction = 0.5_dp*((a-f0)/h+(b-f0)/h)
                weights(i,ib) = weights(i,ib)+length*(1.0_dp-mean_fraction)
                weights(j,ib) = weights(j,ib)+length*mean_fraction
            end do
        end do
        ok = .true.
    end subroutine integration_weights

    function integrate_spectral_bands(spectrum,bands,average) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! CSD (frequency,channel,channel)
        real(dp), intent(in) :: bands(:,:) !! lower and upper frequency edges, (2,nband)
        logical, intent(in), optional :: average !! return mean density instead of integrated band power
        type(multitaper_matrix_result) :: out
        real(dp), allocatable :: weights(:,:)
        integer :: nf, nc, nb, ib, k
        logical :: ok, normalized

        allocate(out%frequencies(0),out%power(0,0,0))
        if (spectrum%status /= spectral_status_success) return
        if (.not. allocated(spectrum%frequencies)) return
        if (.not. allocated(spectrum%power)) return
        nf = size(spectrum%frequencies)
        nc = size(spectrum%power,2)
        if (nf < 2 .or. nc < 2) return
        if (size(spectrum%power,1) /= nf .or. size(spectrum%power,3) /= nc) return
        if (size(bands,1) /= 2) return
        nb = size(bands,2)
        if (nb < 1) return
        if (.not. all(ieee_is_finite(real(spectrum%power,dp)))) return
        if (.not. all(ieee_is_finite(aimag(spectrum%power)))) return
        allocate(weights(nf,nb))
        call integration_weights(spectrum%frequencies,bands,weights,ok)
        if (.not. ok) return
        normalized = .false.
        if (present(average)) normalized = average
        deallocate(out%frequencies,out%power)
        allocate(out%frequencies(nb),out%power(nb,nc,nc))
        out%power = (0.0_dp,0.0_dp)
        do ib = 1, nb
            out%frequencies(ib) = 0.5_dp*bands(1,ib)+0.5_dp*bands(2,ib)
            if (normalized) weights(:,ib) = weights(:,ib)/(bands(2,ib)-bands(1,ib))
            do k = 1, nf
                if (weights(k,ib) <= 0.0_dp) cycle
                out%power(ib,:,:) = out%power(ib,:,:) + &
                    weights(k,ib)*spectrum%power(k,:,:)
            end do
        end do
        if (.not. all(ieee_is_finite(real(out%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(out%power)))) then
            deallocate(out%frequencies,out%power)
            allocate(out%frequencies(0),out%power(0,0,0))
            return
        end if
        ! Cross-frequency observations are not independent tapers: do not
        ! propagate effective_tapers or per-bin predictive shrinkage values.
        out%effective_tapers = 0.0_dp
        out%status = spectral_status_success
    end function integrate_spectral_bands

    function integrate_tf_spectral_bands(spectrum,bands,average) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum !! CSD (f,t,i,j)
        real(dp), intent(in) :: bands(:,:) !! integration bands, (2,nband), in Hz
        logical, intent(in), optional :: average !! normalize by band width when true
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        real(dp), allocatable :: weights(:,:)
        integer :: nf, nt, nc, nb, ib, k
        logical :: ok, normalized

        allocate(out%frequencies(0),out%times(0),out%power(0,0,0,0))
        if (spectrum%status /= spectral_status_success) return
        if (.not. allocated(spectrum%frequencies)) return
        if (.not. allocated(spectrum%times)) return
        if (.not. allocated(spectrum%power)) return
        nf = size(spectrum%frequencies)
        nt = size(spectrum%times)
        nc = size(spectrum%power,3)
        if (nf < 2 .or. nt < 1 .or. nc < 2) return
        if (size(spectrum%power,1) /= nf .or. size(spectrum%power,2) /= nt) return
        if (size(spectrum%power,4) /= nc) return
        if (size(bands,1) /= 2) return
        nb = size(bands,2)
        if (nb < 1) return
        if (.not. all(ieee_is_finite(spectrum%times))) return
        if (.not. all(ieee_is_finite(real(spectrum%power,dp)))) return
        if (.not. all(ieee_is_finite(aimag(spectrum%power)))) return
        allocate(weights(nf,nb))
        call integration_weights(spectrum%frequencies,bands,weights,ok)
        if (.not. ok) return
        normalized = .false.
        if (present(average)) normalized = average
        deallocate(out%frequencies,out%times,out%power)
        allocate(out%frequencies(nb),out%times(nt),out%power(nb,nt,nc,nc))
        out%times = spectrum%times
        out%power = (0.0_dp,0.0_dp)
        do ib = 1, nb
            out%frequencies(ib) = 0.5_dp*bands(1,ib)+0.5_dp*bands(2,ib)
            if (normalized) weights(:,ib) = weights(:,ib)/(bands(2,ib)-bands(1,ib))
            do k = 1, nf
                if (weights(k,ib) <= 0.0_dp) cycle
                out%power(ib,:,:,:) = out%power(ib,:,:,:) + &
                    weights(k,ib)*spectrum%power(k,:,:,:)
            end do
        end do
        if (.not. all(ieee_is_finite(real(out%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(out%power)))) then
            deallocate(out%frequencies,out%times,out%power)
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0,0))
            return
        end if
        out%effective_tapers = 0.0_dp
        out%status = spectral_status_success
    end function integrate_tf_spectral_bands

    function band_multiple_coherence_matrix(spectrum,bands,average,ridge,shrinkage) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! frequency-local CSD matrices
        real(dp), intent(in) :: bands(:,:) !! frequency limits, (2,nband)
        logical, intent(in), optional :: average !! average PSD rather than integrated power
        real(dp), intent(in), optional :: ridge !! relative trace loading, >= 0
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration shrinkage in [0,1]
        type(multiple_coherence_result) :: out
        type(multitaper_matrix_result) :: integrated

        integrated = integrate_spectral_bands(spectrum,bands,average)
        out = multiple_coherence_matrix(integrated,ridge=ridge,shrinkage=shrinkage)
    end function band_multiple_coherence_matrix

    function band_multiple_coherence_tf_matrix(spectrum,bands,average,ridge,shrinkage) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum !! CSD (f,t,i,j)
        real(dp), intent(in) :: bands(:,:) !! frequency limits, (2,nband)
        logical, intent(in), optional :: average !! average PSD rather than integrated power
        real(dp), intent(in), optional :: ridge !! relative trace loading, >= 0
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration shrinkage [0,1]
        type(multiple_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: integrated

        integrated = integrate_tf_spectral_bands(spectrum,bands,average)
        out = multiple_coherence_tf_matrix(integrated,ridge=ridge,shrinkage=shrinkage)
    end function band_multiple_coherence_tf_matrix

    function band_partial_coherence_matrix(spectrum,bands,average,ridge,shrinkage,return_coherency) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! frequency-local CSD matrices
        real(dp), intent(in) :: bands(:,:) !! frequency limits, (2,nband)
        logical, intent(in), optional :: average !! average density rather than integrated power
        real(dp), intent(in), optional :: ridge !! relative trace loading, >= 0
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration shrinkage [0,1]
        logical, intent(in), optional :: return_coherency !! also return phase-bearing coherency
        type(partial_coherence_result) :: out
        type(multitaper_matrix_result) :: integrated

        integrated = integrate_spectral_bands(spectrum,bands,average)
        out = partial_coherence_matrix(integrated,ridge=ridge,shrinkage=shrinkage, &
            return_coherency=return_coherency)
    end function band_partial_coherence_matrix

    function band_partial_coherence_tf_matrix(spectrum,bands,average,ridge,shrinkage,return_coherency) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum !! CSD (f,t,i,j)
        real(dp), intent(in) :: bands(:,:) !! frequency limits, (2,nband)
        logical, intent(in), optional :: average !! average density rather than integrated power
        real(dp), intent(in), optional :: ridge !! relative trace loading, >= 0
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration shrinkage [0,1]
        logical, intent(in), optional :: return_coherency !! also return phase-bearing coherency
        type(partial_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: integrated

        integrated = integrate_tf_spectral_bands(spectrum,bands,average)
        out = partial_coherence_tf_matrix(integrated,ridge=ridge,shrinkage=shrinkage, &
            return_coherency=return_coherency)
    end function band_partial_coherence_tf_matrix

    function multitaper_band_multiple_coherence_real(x,nw,bands,kmax,fs,nfft,detrend, &
            weighting,axis,detrend_mode,average,ridge,shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! observations (sample,channel) unless axis=2
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        real(dp), intent(in) :: bands(:,:) !! requested bands (2,nband) in Hz
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling frequency in Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract mean before tapering
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        logical, intent(in), optional :: average !! output average PSD for each band
        real(dp), intent(in), optional :: ridge !! nonnegative trace ridge
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration shrinkage
        type(multiple_coherence_result) :: out
        type(multitaper_matrix_result) :: spectrum

        spectrum = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_multiple_coherence_matrix(spectrum,bands,average,ridge,shrinkage)
    end function multitaper_band_multiple_coherence_real

    function multitaper_band_multiple_coherence_complex(x,nw,bands,kmax,fs,nfft,detrend, &
            weighting,axis,detrend_mode,average,ridge,shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! observations (sample,channel) unless axis=2
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        real(dp), intent(in) :: bands(:,:) !! requested bands (2,nband) in Hz
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling frequency in Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract mean before tapering
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        logical, intent(in), optional :: average !! output average PSD for each band
        real(dp), intent(in), optional :: ridge !! nonnegative trace ridge
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration shrinkage
        type(multiple_coherence_result) :: out
        type(multitaper_matrix_result) :: spectrum

        spectrum = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_multiple_coherence_matrix(spectrum,bands,average,ridge,shrinkage)
    end function multitaper_band_multiple_coherence_complex

    function multitaper_tf_band_multiple_coherence_real(x,nw,nperseg,bands,noverlap,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! observations (sample,channel) unless axis=2
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! segment sample count
        real(dp), intent(in) :: bands(:,:) !! requested bands (2,nband) in Hz
        integer, intent(in), optional :: noverlap !! number of overlapping samples
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample frequency
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract mean before tapering
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        logical, intent(in), optional :: average !! output mean PSD rather than band power
        real(dp), intent(in), optional :: ridge !! nonnegative trace ridge
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration shrinkage
        type(multiple_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectrum

        spectrum = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_multiple_coherence_tf_matrix(spectrum,bands,average,ridge,shrinkage)
    end function multitaper_tf_band_multiple_coherence_real

    function multitaper_tf_band_multiple_coherence_complex(x,nw,nperseg,bands,noverlap,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! observations (sample,channel) unless axis=2
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! segment sample count
        real(dp), intent(in) :: bands(:,:) !! requested bands (2,nband) in Hz
        integer, intent(in), optional :: noverlap !! number of overlapping samples
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample frequency
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract mean before tapering
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        logical, intent(in), optional :: average !! output mean PSD rather than band power
        real(dp), intent(in), optional :: ridge !! nonnegative trace ridge
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration shrinkage
        type(multiple_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectrum

        spectrum = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_multiple_coherence_tf_matrix(spectrum,bands,average,ridge,shrinkage)
    end function multitaper_tf_band_multiple_coherence_complex

    function multitaper_band_partial_coherence_real(x,nw,bands,kmax,fs,nfft,detrend, &
            weighting,axis,detrend_mode,average,ridge,shrinkage,return_coherency) result(out)
        real(dp), intent(in) :: x(:,:) !! observations (sample,channel) unless axis=2
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        real(dp), intent(in) :: bands(:,:) !! requested frequency bands (2,nband)
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling frequency in Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract sample mean
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! one-based sample axis
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        logical, intent(in), optional :: average !! return average band PSD
        real(dp), intent(in), optional :: ridge !! relative trace loading
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration diagonal shrinkage
        logical, intent(in), optional :: return_coherency !! return complex partial coherency
        type(partial_coherence_result) :: out
        type(multitaper_matrix_result) :: spectrum

        spectrum = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_partial_coherence_matrix(spectrum,bands,average,ridge,shrinkage,return_coherency)
    end function multitaper_band_partial_coherence_real

    function multitaper_tf_band_partial_coherence_real(x,nw,nperseg,bands,noverlap,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage,return_coherency) result(out)
        real(dp), intent(in) :: x(:,:) !! observations (sample,channel) unless axis=2
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! segment length in samples
        real(dp), intent(in) :: bands(:,:) !! requested frequency bands (2,nband)
        integer, intent(in), optional :: noverlap !! overlapping sample count
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! positive sampling rate in Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract segment mean
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! one-based sample axis
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        logical, intent(in), optional :: average !! return average band PSD
        real(dp), intent(in), optional :: ridge !! relative trace loading
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration diagonal shrinkage
        logical, intent(in), optional :: return_coherency !! return complex partial coherency
        type(partial_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectrum

        spectrum = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_partial_coherence_tf_matrix(spectrum,bands,average,ridge,shrinkage,return_coherency)
    end function multitaper_tf_band_partial_coherence_real

    function multitaper_band_partial_coherence_complex(x,nw,bands,kmax,fs,nfft,detrend, &
            weighting,axis,detrend_mode,average,ridge,shrinkage,return_coherency) result(out)
        complex(dp), intent(in) :: x(:,:) !! observations (sample,channel) unless axis=2
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        real(dp), intent(in) :: bands(:,:) !! requested frequency bands (2,nband)
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling frequency in Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract sample mean
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! one-based sample axis
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        logical, intent(in), optional :: average !! return average band PSD
        real(dp), intent(in), optional :: ridge !! relative trace loading
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration diagonal shrinkage
        logical, intent(in), optional :: return_coherency !! return complex partial coherency
        type(partial_coherence_result) :: out
        type(multitaper_matrix_result) :: spectrum

        spectrum = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_partial_coherence_matrix(spectrum,bands,average,ridge,shrinkage,return_coherency)
    end function multitaper_band_partial_coherence_complex

    function multitaper_tf_band_partial_coherence_complex(x,nw,nperseg,bands,noverlap,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage,return_coherency) result(out)
        complex(dp), intent(in) :: x(:,:) !! observations (sample,channel) unless axis=2
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! segment length in samples
        real(dp), intent(in) :: bands(:,:) !! requested frequency bands (2,nband)
        integer, intent(in), optional :: noverlap !! overlapping sample count
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! positive sampling rate in Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract segment mean
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! one-based sample axis
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        logical, intent(in), optional :: average !! return average band PSD
        real(dp), intent(in), optional :: ridge !! relative trace loading
        real(dp), intent(in), optional :: shrinkage !! explicit post-integration diagonal shrinkage
        logical, intent(in), optional :: return_coherency !! return complex partial coherency
        type(partial_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectrum

        spectrum = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_partial_coherence_tf_matrix(spectrum,bands,average,ridge,shrinkage,return_coherency)
    end function multitaper_tf_band_partial_coherence_complex
end module scifort_signal_band_coherence
