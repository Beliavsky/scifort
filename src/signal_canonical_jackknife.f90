! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Window-delete jackknife for conditional canonical spectral coherence.
! The resampling units are nonoverlapping time windows, or contiguous
! equal-length blocks thereof, never neighboring frequency bins or tapers.
! Confidence limits are pointwise, heuristic normal-approximation intervals
! on the Fisher-atanh(sqrt(coherence)) scale, NOT calibrated significance tests.
! This is a SciFort extension, not a scipy.signal function.
module scifort_signal_canonical_jackknife
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix_spectrogram
    use scifort_signal_band_coherence, only : integrate_tf_spectral_bands
    use scifort_signal_canonical_coherence, only : canonical_coherence_matrix, canonical_coherence_result
    use scifort_signal_partial_coherence, only : partial_coherence_success, &
        partial_coherence_invalid_input, partial_coherence_singular
    implicit none
    private

    real(dp), parameter :: normal_975 = 1.959963984540054_dp
    real(dp), parameter :: fisher_ceiling = 1.0_dp - 1.0e-12_dp

    type, public :: canonical_window_jackknife_result
        real(dp), allocatable :: frequencies(:) !! frequencies or band centers
        real(dp), allocatable :: squared_modes(:,:) !! full-mean canonical roots (frequency,rank)
        real(dp), allocatable :: stderr(:,:) !! delete-one-window jackknife SE of squared roots
        real(dp), allocatable :: stderr_fisher(:,:) !! SE on atanh(sqrt(root)) scale
        real(dp), allocatable :: lower95(:,:) !! heuristic 95% pointwise Fisher-scale lower limits
        real(dp), allocatable :: upper95(:,:) !! heuristic 95% pointwise Fisher-scale upper limits
        real(dp), allocatable :: bias(:,:) !! estimated jackknife bias of squared roots
        logical, allocatable :: valid(:,:) !! inference valid for each isolated interior root
        real(dp), allocatable :: coherence(:) !! full-mean aggregate group coherence
        real(dp), allocatable :: coherence_stderr(:) !! jackknife SE of group coherence
        logical, allocatable :: coherence_valid(:) !! all leave-one-out aggregate fits valid
        integer :: n_windows = 0 !! total number of input nonoverlapping windows
        integer :: n_blocks = 0 !! number of jackknife deletion groups (>=4 on success)
        integer :: block_length = 0 !! number of consecutive windows deleted per group
        integer :: status = partial_coherence_invalid_input !! invalid, success or singular
    end type canonical_window_jackknife_result

    interface multitaper_canonical_window_jackknife
        module procedure multitaper_window_jackknife_real
        module procedure multitaper_window_jackknife_complex
    end interface

    interface multitaper_canonical_block_jackknife
        module procedure multitaper_block_jackknife_real
        module procedure multitaper_block_jackknife_complex
    end interface

    public :: canonical_window_jackknife, multitaper_canonical_window_jackknife
    public :: canonical_block_jackknife, multitaper_canonical_block_jackknife

contains

    pure real(dp) function to_fisher(squared) result(value)
        real(dp), intent(in) :: squared !! squared coherence in [0,1]
        value = atanh(sqrt(max(0.0_dp,min(fisher_ceiling,squared))))
    end function to_fisher

    pure logical function isolated_mode(modes,rank,tolerance) result(ok)
        real(dp), intent(in) :: modes(:) !! descending canonical squared roots
        integer, intent(in) :: rank !! mode index, 1-based
        real(dp), intent(in) :: tolerance !! minimum separation from neighbors/boundaries

        ok = .false.
        if (.not. all(ieee_is_finite(modes))) return
        if (modes(rank) <= tolerance .or. modes(rank) >= 1.0_dp-tolerance) return
        if (rank > 1) then
            if (modes(rank-1)-modes(rank) <= tolerance) return
        end if
        if (rank < size(modes)) then
            if (modes(rank)-modes(rank+1) <= tolerance) return
        end if
        ok = .true.
    end function isolated_mode

    function canonical_window_jackknife(spectra,group_a,group_b,conditioned_on,ridge,shrinkage, &
            bands,average,gap_tolerance,block_length) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: &
            spectra !! complete CSDs (frequency,window,channel,channel)
        integer, intent(in) :: group_a(:) !! first nonempty group of 1-based channels
        integer, intent(in) :: group_b(:) !! second nonempty group of disjoint channels
        integer, intent(in), optional :: conditioned_on(:) !! channels to condition on
        real(dp), intent(in), optional :: ridge !! trace-relative diagonal loading, >=0
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal target coefficient, [0,1]
        real(dp), intent(in), optional :: bands(:,:) !! optional (2,nband) integration bounds
        logical, intent(in), optional :: average !! average rather than integrate band CSD
        real(dp), intent(in), optional :: gap_tolerance !! ignore boundary/tied roots (default 1e-6)
        integer, intent(in), optional :: block_length !! consecutive deleted windows per block (default 1)
        type(canonical_window_jackknife_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: work
        type(multitaper_matrix_result) :: mean_spectrum
        type(canonical_coherence_result) :: point, leave
        real(dp), allocatable :: leave_roots(:,:,:), leave_groups(:)
        complex(dp), allocatable :: totals(:,:,:)
        real(dp) :: nan, tol, avg, se, mean_z, se_z, center_z
        integer :: nfreq, nwin, nchan, nm, nblocks, bl, f, mode, k, first, last
        logical :: regular, all_good

        if (spectra%status /= partial_coherence_success) return
        if (.not. allocated(spectra%power) .or. .not. allocated(spectra%times) .or. &
            .not. allocated(spectra%frequencies)) return
        nwin = size(spectra%times)
        nfreq = size(spectra%frequencies)
        nchan = size(spectra%power,3)
        bl = 1
        if (present(block_length)) bl = block_length
        if (bl < 1) return
        if (mod(nwin,bl) /= 0) return
        nblocks = nwin/bl
        if (nblocks < 4 .or. nfreq < 1 .or. nchan < 2) return
        if (size(spectra%power,1) /= nfreq .or. size(spectra%power,2) /= nwin .or. &
            size(spectra%power,4) /= nchan) return
        if (.not. all(ieee_is_finite(spectra%times))) return
        if (any(spectra%times(2:) <= spectra%times(:nwin-1))) return
        if (.not. all(ieee_is_finite(spectra%frequencies))) return
        if (.not. all(ieee_is_finite(real(spectra%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(spectra%power)))) return
        tol = 1.0e-6_dp
        if (present(gap_tolerance)) tol = gap_tolerance
        if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp .or. tol >= 0.5_dp) return

        if (present(bands)) then
            work = integrate_tf_spectral_bands(spectra,bands,average)
        else
            work = spectra
        end if
        if (work%status /= partial_coherence_success) return
        nfreq = size(work%frequencies)
        allocate(mean_spectrum%frequencies(nfreq),mean_spectrum%power(nfreq,nchan,nchan), &
            totals(nfreq,nchan,nchan))
        totals = sum(work%power,dim=2)
        mean_spectrum%frequencies = work%frequencies
        mean_spectrum%power = totals/real(nwin,dp)
        mean_spectrum%status = partial_coherence_success
        ! Each delete-one fit uses the exact same fixed regularizer; auto/CV
        ! shrinkage across adaptive tapers would require the raw tapered FFTs.
        point = canonical_coherence_matrix(mean_spectrum,group_a,group_b, &
            conditioned_on,ridge=ridge,shrinkage=shrinkage)
        if (point%status == partial_coherence_invalid_input) return
        nm = size(point%squared_modes,2)
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        allocate(out%frequencies(nfreq),out%squared_modes(nfreq,nm),out%stderr(nfreq,nm), &
            out%stderr_fisher(nfreq,nm),out%lower95(nfreq,nm),out%upper95(nfreq,nm), &
            out%bias(nfreq,nm),out%valid(nfreq,nm),out%coherence(nfreq), &
            out%coherence_stderr(nfreq),out%coherence_valid(nfreq), &
            leave_roots(nfreq,nm,nblocks),leave_groups(nblocks))
        out%frequencies = point%frequencies
        out%squared_modes = point%squared_modes
        out%coherence = point%coherence
        out%stderr = nan
        out%stderr_fisher = nan
        out%lower95 = nan
        out%upper95 = nan
        out%bias = nan
        out%valid = .false.
        out%coherence_stderr = nan
        out%coherence_valid = .false.
        out%n_windows = nwin
        out%n_blocks = nblocks
        out%block_length = bl
        out%status = point%status
        leave_roots = nan
        leave_groups = nan
        ! Delete disjoint, contiguous, equal-sized blocks of complete windows.
        ! The standard equal-group jackknife multiplier uses nblocks, not nwin.
        do k = 1, nblocks
            first = (k-1)*bl+1
            last = k*bl
            mean_spectrum%power = (totals-sum(work%power(:,first:last,:,:),dim=2))/real(nwin-bl,dp)
            leave = canonical_coherence_matrix(mean_spectrum,group_a,group_b, &
                conditioned_on,ridge=ridge,shrinkage=shrinkage)
            if (leave%status == partial_coherence_invalid_input) then
                out%status = partial_coherence_invalid_input
                return
            end if
            do f = 1, nfreq
                if (.not. leave%valid(f)) cycle
                leave_roots(f,:,k) = leave%squared_modes(f,:)
            end do
        end do
        do f = 1, nfreq
            if (.not. point%valid(f)) cycle
            ! All group coherence statistics are functions of the canonical roots.
            all_good = all(ieee_is_finite(leave_roots(f,:,:)))
            if (.not. all_good) cycle
            do k = 1, nblocks
                leave_groups(k) = 1.0_dp-product(1.0_dp-leave_roots(f,:,k))
            end do
            avg = sum(leave_groups)/real(nblocks,dp)
            se = sqrt(real(nblocks-1,dp)/real(nblocks,dp)* &
                sum((leave_groups-avg)**2))
            out%coherence_stderr(f) = se
            out%coherence_valid(f) = .true.
            do mode = 1, nm
                regular = isolated_mode(point%squared_modes(f,:),mode,tol)
                do k = 1, nblocks
                    if (.not. isolated_mode(leave_roots(f,:,k),mode,tol)) regular = .false.
                end do
                if (.not. regular) cycle
                avg = sum(leave_roots(f,mode,:))/real(nblocks,dp)
                out%bias(f,mode) = real(nblocks-1,dp)*(avg-point%squared_modes(f,mode))
                out%stderr(f,mode) = sqrt(real(nblocks-1,dp)/real(nblocks,dp)* &
                    sum((leave_roots(f,mode,:)-avg)**2))
                mean_z = 0.0_dp
                do k = 1, nblocks
                    mean_z = mean_z+to_fisher(leave_roots(f,mode,k))
                end do
                mean_z = mean_z/real(nblocks,dp)
                se_z = 0.0_dp
                do k = 1, nblocks
                    se_z = se_z+(to_fisher(leave_roots(f,mode,k))-mean_z)**2
                end do
                se_z = sqrt(real(nblocks-1,dp)/real(nblocks,dp)*se_z)
                center_z = to_fisher(point%squared_modes(f,mode))
                out%stderr_fisher(f,mode) = se_z
                out%lower95(f,mode) = tanh(max(0.0_dp,center_z-normal_975*se_z))**2
                out%upper95(f,mode) = tanh(center_z+normal_975*se_z)**2
                out%valid(f,mode) = .true.
            end do
        end do
        if (.not. all(out%coherence_valid)) out%status = partial_coherence_singular
    end function canonical_window_jackknife

    function multitaper_window_jackknife_real(x,nw,nperseg,group_a,group_b,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            bands,average,gap_tolerance,block_length) result(out)
        real(dp), intent(in) :: x(:,:) !! samples by channel (or axis-selected layout)
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples per independent, nonoverlapping window
        integer, intent(in) :: group_a(:) !! group A channel indices
        integer, intent(in) :: group_b(:) !! group B channel indices
        integer, intent(in), optional :: conditioned_on(:) !! conditional channels
        integer, intent(in), optional :: kmax !! DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling rate, Hz
        integer, intent(in), optional :: nfft !! FFT size
        logical, intent(in), optional :: detrend !! remove segment mean
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! time axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant or linear
        real(dp), intent(in), optional :: ridge !! fractional diagonal loading
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage
        real(dp), intent(in), optional :: bands(:,:) !! optional frequency integration bands
        logical, intent(in), optional :: average !! average CSD in each band
        real(dp), intent(in), optional :: gap_tolerance !! near-degenerate mode threshold
        integer, intent(in), optional :: block_length !! consecutive deleted windows per jackknife block
        type(canonical_window_jackknife_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectral_windows

        spectral_windows = multitaper_spectral_matrix_spectrogram(x,nw,nperseg, &
            noverlap=0,kmax=kmax,fs=fs,nfft=nfft,detrend=detrend, &
            weighting=weighting,axis=axis,detrend_mode=detrend_mode)
        out = canonical_window_jackknife(spectral_windows,group_a,group_b, &
            conditioned_on,ridge,shrinkage,bands,average,gap_tolerance,block_length)
    end function multitaper_window_jackknife_real

    function multitaper_window_jackknife_complex(x,nw,nperseg,group_a,group_b,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            bands,average,gap_tolerance,block_length) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex time-series samples
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! independent, nonoverlapping sample window size
        integer, intent(in) :: group_a(:) !! first group
        integer, intent(in) :: group_b(:) !! second group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: kmax !! DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling rate
        integer, intent(in), optional :: nfft !! FFT size
        logical, intent(in), optional :: detrend !! remove segment mean
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! time axis
        character(*), intent(in), optional :: detrend_mode !! detrending mode
        real(dp), intent(in), optional :: ridge !! fractional trace loading
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage
        real(dp), intent(in), optional :: bands(:,:) !! integration bounds
        logical, intent(in), optional :: average !! mean within band instead of integral
        real(dp), intent(in), optional :: gap_tolerance !! reject nearly tied roots
        integer, intent(in), optional :: block_length !! consecutive deleted windows per jackknife block
        type(canonical_window_jackknife_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectral_windows

        spectral_windows = multitaper_spectral_matrix_spectrogram(x,nw,nperseg, &
            noverlap=0,kmax=kmax,fs=fs,nfft=nfft,detrend=detrend, &
            weighting=weighting,axis=axis,detrend_mode=detrend_mode)
        out = canonical_window_jackknife(spectral_windows,group_a,group_b, &
            conditioned_on,ridge,shrinkage,bands,average,gap_tolerance,block_length)
    end function multitaper_window_jackknife_complex

    function canonical_block_jackknife(spectra,group_a,group_b,block_length,conditioned_on, &
            ridge,shrinkage,bands,average,gap_tolerance) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectra !! input window spectra
        integer, intent(in) :: group_a(:) !! first group of channel indices
        integer, intent(in) :: group_b(:) !! second group of channel indices
        integer, intent(in) :: block_length !! equally sized, consecutive window blocks, >=1
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channel indices
        real(dp), intent(in), optional :: ridge !! trace-relative diagonal loading
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage
        real(dp), intent(in), optional :: bands(:,:) !! frequency integration bounds (2,nband)
        logical, intent(in), optional :: average !! mean, rather than integral, for each band
        real(dp), intent(in), optional :: gap_tolerance !! mode separation threshold
        type(canonical_window_jackknife_result) :: out

        out = canonical_window_jackknife(spectra,group_a,group_b,conditioned_on, &
            ridge,shrinkage,bands,average,gap_tolerance,block_length)
    end function canonical_block_jackknife

    function multitaper_block_jackknife_real(x,nw,nperseg,group_a,group_b,block_length,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            bands,average,gap_tolerance) result(out)
        real(dp), intent(in) :: x(:,:) !! real signals, samples by channel or axis-selected layout
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples in each nonoverlapping spectral window
        integer, intent(in) :: group_a(:) !! group A channel indices
        integer, intent(in) :: group_b(:) !! group B channel indices
        integer, intent(in) :: block_length !! consecutive whole spectral windows per deletion group
        integer, intent(in), optional :: conditioned_on(:) !! conditional group indices
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sample frequency
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove constant from each window
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! time dimension, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        real(dp), intent(in), optional :: ridge !! trace-relative diagonal loading
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        real(dp), intent(in), optional :: bands(:,:) !! optional band limits (2,nband)
        logical, intent(in), optional :: average !! use mean CSD instead of integral within a band
        real(dp), intent(in), optional :: gap_tolerance !! reject nearly tied/boundary roots
        type(canonical_window_jackknife_result) :: out

        out = multitaper_window_jackknife_real(x,nw,nperseg,group_a,group_b,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            bands,average,gap_tolerance,block_length)
    end function multitaper_block_jackknife_real

    function multitaper_block_jackknife_complex(x,nw,nperseg,group_a,group_b,block_length,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            bands,average,gap_tolerance) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex signals, samples by channel or axis-selected layout
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples in each nonoverlapping spectral window
        integer, intent(in) :: group_a(:) !! group A channel indices
        integer, intent(in) :: group_b(:) !! group B channel indices
        integer, intent(in) :: block_length !! consecutive whole spectral windows per deletion group
        integer, intent(in), optional :: conditioned_on(:) !! conditional group indices
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sample frequency
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove constant from each window
        character(*), intent(in), optional :: weighting !! equal, eigen, adaptive
        integer, intent(in), optional :: axis !! time dimension, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, linear
        real(dp), intent(in), optional :: ridge !! trace-relative diagonal loading
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        real(dp), intent(in), optional :: bands(:,:) !! optional band limits (2,nband)
        logical, intent(in), optional :: average !! use mean CSD instead of integral within a band
        real(dp), intent(in), optional :: gap_tolerance !! reject nearly tied/boundary roots
        type(canonical_window_jackknife_result) :: out

        out = multitaper_window_jackknife_complex(x,nw,nperseg,group_a,group_b,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            bands,average,gap_tolerance,block_length)
    end function multitaper_block_jackknife_complex

end module scifort_signal_canonical_jackknife
