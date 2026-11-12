! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Exploratory, autocovariance-pilot selection of the moving-block length for
! canonical spectral coherence. This is a SciFort extension, not scipy.signal.
! The plug-in scaling is a heuristic for dependent window-level CSD estimates;
! it is not an optimal-length or calibrated coverage guarantee.
module scifort_signal_canonical_block_selection
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_multitaper, only : multitaper_matrix_result
    use scifort_signal_band_coherence, only : integrate_tf_spectral_bands
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, canonical_coherence_matrix
    use scifort_signal_partial_coherence, only : partial_coherence_success, partial_coherence_invalid_input
    implicit none
    private

    type, public :: canonical_block_length_result
        integer :: block_length = 0 !! proposed moving-block length in spectral windows
        integer :: n_windows = 0 !! complete spectral windows supplied
        integer :: max_lag = 0 !! maximum autocovariance lag evaluated
        integer :: lags_used = 0 !! successive detectable autocorrelations in the pilot
        integer :: maximum_block_length = 0 !! upper limit applied to selected length
        real(dp) :: first_lag_correlation = 0.0_dp !! normalized lag-one Frobenius autocovariance
        real(dp) :: long_run_variance_ratio = 1.0_dp !! tapered long-run variance / lag-zero variance
        real(dp) :: weighted_lag_ratio = 0.0_dp !! tapered first lag moment / long-run variance
        real(dp), allocatable :: autocorrelations(:) !! pilot normalized correlations, lags 1:max_lag
        integer :: mode_index = 0 !! 0: group coherence; positive: chosen squared canonical mode
        integer :: frequency_index = 0 !! 0: mean over frequencies; positive: one frequency/bin
        real(dp) :: full_statistic = 0.0_dp !! full-sample target statistic for influence pilot
        real(dp), allocatable :: influence(:) !! delete-one influence scores, window order (optional)
        integer :: status = partial_coherence_invalid_input !! success or invalid input
    end type canonical_block_length_result

    public :: select_canonical_block_length, select_canonical_influence_block_length

contains

    function select_canonical_block_length(spectra,max_lag,max_block_length) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectra !! Hermitian CSD (freq,time,channel,channel)
        integer, intent(in), optional :: max_lag !! pilot lag cap, from 1 through n_windows-2
        integer, intent(in), optional :: max_block_length !! permitted length cap, from 1 through n_windows-1
        type(canonical_block_length_result) :: out
        complex(dp), allocatable :: centered(:,:,:,:), mean_matrix(:,:,:)
        real(dp) :: scale, gamma0, gamma, cutoff, long_run, moment, weight, ratio, pilot_length
        integer :: nf, nw, nc, k, lag_cap, length_cap, nused, t

        if (spectra%status /= partial_coherence_success) return
        if (.not. allocated(spectra%power) .or. .not. allocated(spectra%times) .or. &
            .not. allocated(spectra%frequencies)) return
        nf = size(spectra%frequencies)
        nw = size(spectra%times)
        nc = size(spectra%power,3)
        if (nf < 1 .or. nw < 4 .or. nc < 2) return
        if (size(spectra%power,1) /= nf .or. size(spectra%power,2) /= nw .or. &
            size(spectra%power,4) /= nc) return
        if (.not. all(ieee_is_finite(spectra%times)) .or. &
            .not. all(ieee_is_finite(spectra%frequencies))) return
        if (any(spectra%times(2:) <= spectra%times(:nw-1))) return
        if (.not. all(ieee_is_finite(real(spectra%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(spectra%power)))) return

        lag_cap = min(nw-2,max(1,int(sqrt(real(nw,dp)))))
        if (present(max_lag)) lag_cap = max_lag
        if (lag_cap < 1 .or. lag_cap > nw-2) return
        length_cap = max(1,min(nw-1,nw/3))
        if (present(max_block_length)) length_cap = max_block_length
        if (length_cap < 1 .or. length_cap >= nw) return

        out%n_windows = nw
        out%max_lag = lag_cap
        out%maximum_block_length = length_cap
        out%block_length = 1
        out%status = partial_coherence_success
        allocate(out%autocorrelations(lag_cap))
        out%autocorrelations = 0.0_dp

        ! A common overall scale is removed before taking products, so large
        ! spectral units do not overflow the lag-covariance calculation.
        allocate(mean_matrix(nf,nc,nc),centered(nf,nw,nc,nc))
        scale = maxval(abs(spectra%power))
        if (scale == 0.0_dp) return
        centered = spectra%power/scale
        mean_matrix = sum(centered,dim=2)/real(nw,dp)
        do t = 1, nw
            centered(:,t,:,:) = centered(:,t,:,:)-mean_matrix
        end do
        scale = maxval(abs(centered))
        if (scale == 0.0_dp) return
        if (.not. ieee_is_finite(scale)) then
            out%status = partial_coherence_invalid_input
            out%block_length = 0
            return
        end if
        centered = centered/scale
        gamma0 = real(sum(conjg(centered)*centered),dp)/real(nw,dp)
        if (gamma0 <= 0.0_dp) return

        ! Treat each complex Hermitian matrix as a real vector with Frobenius
        ! inner product. This aggregates temporal dependence without assuming
        ! that different frequencies or channel entries are independent.
        cutoff = 2.0_dp/sqrt(real(nw,dp))
        nused = 0
        do k = 1, lag_cap
            gamma = real(sum(conjg(centered(:,1:nw-k,:,:))* &
                centered(:,1+k:nw,:,:)),dp)/real(nw-k,dp)
            out%autocorrelations(k) = gamma/gamma0
            if (k == 1) out%first_lag_correlation = out%autocorrelations(k)
            if (abs(out%autocorrelations(k)) <= cutoff) exit
            nused = k
        end do
        out%lags_used = nused
        if (nused == 0) return

        ! Signed autocovariances with a Bartlett taper estimate the long-run
        ! variance and its first weighted lag moment. The n^(1/3) scaling is
        ! appropriate for variance-oriented block-size plug-in heuristics,
        ! but the pilot cutoff and constants are deliberately exploratory.
        long_run = 1.0_dp
        moment = 0.0_dp
        do k = 1, nused
            weight = 1.0_dp-real(k,dp)/real(nused+1,dp)
            long_run = long_run+2.0_dp*weight*out%autocorrelations(k)
            moment = moment+2.0_dp*real(k,dp)*weight*out%autocorrelations(k)
        end do
        if (.not. ieee_is_finite(long_run) .or. .not. ieee_is_finite(moment)) then
            out%status = partial_coherence_invalid_input
            out%block_length = 0
            return
        end if
        if (long_run <= sqrt(epsilon(1.0_dp))) return
        out%long_run_variance_ratio = long_run
        ratio = moment/long_run
        out%weighted_lag_ratio = ratio
        pilot_length = (1.5_dp*real(nw,dp)*ratio**2)**(1.0_dp/3.0_dp)
        if (.not. ieee_is_finite(pilot_length)) then
            out%status = partial_coherence_invalid_input
            out%block_length = 0
            return
        end if
        out%block_length = max(1,min(length_cap,nint(min(real(length_cap,dp),pilot_length))))
    end function select_canonical_block_length

    ! Numerically evaluate delete-one-window influence scores for the *requested*
    ! canonical statistic. The leave-one estimates use the same fixed ridge and
    ! shrinkage as the final bootstrap, with no data-dependent re-tuning.
    ! This is an exploratory block-length pilot, not a formal influence-function
    ! theorem. It never treats adjacent frequencies as independent observations.
    function select_canonical_influence_block_length(spectra,group_a,group_b, &
            conditioned_on,ridge,shrinkage,bands,average,mode_index,frequency_index, &
            max_lag,max_block_length) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectra !! window CSDs (freq,window,channel,channel)
        integer, intent(in) :: group_a(:) !! one-based nonempty first group
        integer, intent(in) :: group_b(:) !! one-based nonempty disjoint second group
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels
        real(dp), intent(in), optional :: ridge !! fixed nonnegative ridge loading
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage in [0,1]
        real(dp), intent(in), optional :: bands(:,:) !! optional spectral bands (2,nbands)
        logical, intent(in), optional :: average !! average instead of integrate band matrices
        integer, intent(in), optional :: mode_index !! 0 for group coherence; 1..nmode for squared canonical root
        integer, intent(in), optional :: frequency_index !! 0 for all-frequency average; 1..nf for one bin
        integer, intent(in), optional :: max_lag !! pilot maximum autocorrelation lag
        integer, intent(in), optional :: max_block_length !! upper cap on selected block length
        type(canonical_block_length_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: work, proxy
        type(multitaper_matrix_result) :: mean_spectrum
        type(canonical_coherence_result) :: full, omitted
        complex(dp), allocatable :: total(:,:,:)
        real(dp), allocatable :: influence(:)
        real(dp) :: full_value, leave_value
        integer :: nw, nf, nc, t, f0, f1, mode, selected_frequency

        if (spectra%status /= partial_coherence_success) return
        if (.not. allocated(spectra%power) .or. .not. allocated(spectra%times) .or. &
            .not. allocated(spectra%frequencies)) return
        nw = size(spectra%times)
        nf = size(spectra%frequencies)
        nc = size(spectra%power,3)
        if (nw < 4 .or. nf < 1 .or. nc < 2) return
        if (size(spectra%power,1) /= nf .or. size(spectra%power,2) /= nw .or. &
            size(spectra%power,4) /= nc) return
        if (.not. all(ieee_is_finite(spectra%times)) .or. &
            .not. all(ieee_is_finite(spectra%frequencies))) return
        if (any(spectra%times(2:) <= spectra%times(:nw-1))) return
        if (.not. all(ieee_is_finite(real(spectra%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(spectra%power)))) return

        mode = 0
        if (present(mode_index)) mode = mode_index
        if (mode < 0 .or. mode > min(size(group_a),size(group_b))) return
        selected_frequency = 0
        if (present(frequency_index)) selected_frequency = frequency_index

        if (present(bands)) then
            work = integrate_tf_spectral_bands(spectra,bands,average)
        else
            work = spectra
        end if
        if (work%status /= partial_coherence_success) return
        nf = size(work%frequencies)
        if (selected_frequency < 0 .or. selected_frequency > nf) return
        f0 = 1
        f1 = nf
        if (selected_frequency > 0) then
            f0 = selected_frequency
            f1 = selected_frequency
        end if

        allocate(total(nf,nc,nc),mean_spectrum%frequencies(nf), &
            mean_spectrum%power(nf,nc,nc),influence(nw))
        mean_spectrum%frequencies = work%frequencies
        mean_spectrum%status = partial_coherence_success
        total = sum(work%power,dim=2)
        mean_spectrum%power = total/real(nw,dp)
        full = canonical_coherence_matrix(mean_spectrum,group_a,group_b, &
            conditioned_on,ridge=ridge,shrinkage=shrinkage)
        if (full%status == partial_coherence_invalid_input) return
        if (.not. allocated(full%valid)) return
        if (.not. all(full%valid(f0:f1))) return
        if (mode > size(full%squared_modes,2)) return
        if (mode > 0) then
            if (.not. isolated_target(full%squared_modes(f0:f1,:),mode)) return
        end if
        full_value = target_value(full,mode,f0,f1)
        if (.not. ieee_is_finite(full_value)) return

        do t = 1, nw
            mean_spectrum%power = (total-work%power(:,t,:,:))/real(nw-1,dp)
            omitted = canonical_coherence_matrix(mean_spectrum,group_a,group_b, &
                conditioned_on,ridge=ridge,shrinkage=shrinkage)
            if (omitted%status == partial_coherence_invalid_input) return
            if (.not. allocated(omitted%valid)) return
            if (.not. all(omitted%valid(f0:f1))) return
            if (mode > 0) then
                if (.not. isolated_target(omitted%squared_modes(f0:f1,:),mode)) return
            end if
            leave_value = target_value(omitted,mode,f0,f1)
            if (.not. ieee_is_finite(leave_value)) return
            influence(t) = real(nw-1,dp)*(full_value-leave_value)
        end do
        if (.not. all(ieee_is_finite(influence))) return

        ! The established Frobenius pilot, applied to a single real series,
        ! yields *exactly* the scalar autocovariances of these influences.
        ! The dummy second channel is zero; no spectral PD check is required
        ! because this proxy is used strictly to select the block length.
        allocate(proxy%power(1,nw,2,2),proxy%times(nw),proxy%frequencies(1))
        proxy%status = partial_coherence_success
        proxy%power = (0.0_dp,0.0_dp)
        proxy%power(1,:,1,1) = cmplx(influence,0.0_dp,dp)
        proxy%times = spectra%times
        proxy%frequencies = [0.0_dp]
        out = select_canonical_block_length(proxy,max_lag,max_block_length)
        if (out%status /= partial_coherence_success) return
        out%mode_index = mode
        out%frequency_index = selected_frequency
        out%full_statistic = full_value
        call move_alloc(influence,out%influence)
    end function select_canonical_influence_block_length

    pure real(dp) function target_value(fit,mode,f0,f1) result(value)
        type(canonical_coherence_result), intent(in) :: fit !! valid canonical fit for selected bins
        integer, intent(in) :: mode !! zero for group coherence; otherwise canonical mode
        integer, intent(in) :: f0 !! first selected frequency
        integer, intent(in) :: f1 !! last selected frequency

        if (mode == 0) then
            value = sum(fit%coherence(f0:f1))/real(f1-f0+1,dp)
        else
            value = sum(fit%squared_modes(f0:f1,mode))/real(f1-f0+1,dp)
        end if
    end function target_value

    pure logical function isolated_target(modes,mode) result(ok)
        real(dp), intent(in) :: modes(:,:) !! squared roots for selected bins, descending
        integer, intent(in) :: mode !! one-based target root
        integer :: f
        real(dp), parameter :: gap = 1.0e-6_dp

        ok = .false.
        do f = 1, size(modes,1)
            if (.not. all(ieee_is_finite(modes(f,:)))) return
            if (modes(f,mode) <= gap .or. modes(f,mode) >= 1.0_dp-gap) return
            if (mode > 1) then
                if (modes(f,mode-1)-modes(f,mode) <= gap) return
            end if
            if (mode < size(modes,2)) then
                if (modes(f,mode)-modes(f,mode+1) <= gap) return
            end if
        end do
        ok = .true.
    end function isolated_target

end module scifort_signal_canonical_block_selection
