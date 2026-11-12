! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Moving-block bootstrap for canonical spectral coherence. Complete,
! nonoverlapping spectral-estimation windows are the resampling units.
! This SciFort extension is not part of scipy.signal. Its bootstrap
! intervals are exploratory pointwise summaries, not calibrated tests.
module scifort_signal_canonical_bootstrap
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix_spectrogram
    use scifort_signal_band_coherence, only : integrate_tf_spectral_bands
    use scifort_signal_canonical_coherence, only : canonical_coherence_matrix, canonical_coherence_result
    use scifort_signal_canonical_block_selection, only : canonical_block_length_result, &
        select_canonical_block_length, select_canonical_influence_block_length
    use scifort_signal_partial_coherence, only : partial_coherence_success, &
        partial_coherence_invalid_input, partial_coherence_singular
    implicit none
    private

    integer(int64), parameter :: rng_modulus = 2147483647_int64
    integer(int64), parameter :: rng_multiplier = 16807_int64
    integer, parameter :: default_resamples = 399
    integer, parameter :: default_seed = 12345
    real(dp), parameter :: normal_975 = 1.959963984540054_dp

    type, public :: canonical_moving_block_bootstrap_result
        real(dp), allocatable :: frequencies(:) !! frequencies or integrated band centers
        real(dp), allocatable :: squared_modes(:,:) !! full-sample roots (frequency,mode)
        real(dp), allocatable :: stderr(:,:) !! bootstrap sample standard deviations of squared roots
        real(dp), allocatable :: bias(:,:) !! bootstrap mean minus original squared roots
        real(dp), allocatable :: studentizing_stderr(:,:) !! original delete-one-block standard errors; allocated for studentized
        real(dp), allocatable :: bca_acceleration(:,:) !! block-jackknife acceleration, allocated only for bca
        real(dp), allocatable :: lower95(:,:) !! 95% bounds using selected interval method
        real(dp), allocatable :: upper95(:,:) !! selected-method upper 95% bounds
        logical, allocatable :: valid(:,:) !! all refits successful, roots isolated and interior
        real(dp), allocatable :: coherence(:) !! aggregate full-sample group coherence
        real(dp), allocatable :: coherence_stderr(:) !! bootstrap aggregate standard deviation
        real(dp), allocatable :: coherence_bias(:) !! bootstrap mean minus full-sample aggregate
        real(dp), allocatable :: coherence_studentizing_stderr(:) !! original block-jackknife standard errors
        real(dp), allocatable :: coherence_bca_acceleration(:) !! aggregate block acceleration, allocated only for bca
        real(dp), allocatable :: coherence_lower95(:) !! selected-method aggregate lower bound
        real(dp), allocatable :: coherence_upper95(:) !! selected-method aggregate upper bound
        logical, allocatable :: coherence_valid(:) !! all bootstrap aggregate refits succeeded
        integer :: n_windows = 0 !! number of complete spectral-estimation windows
        integer :: n_resamples = 0 !! number of independently drawn bootstrap replicates
        integer :: block_length = 0 !! consecutive windows in each drawn block
        integer :: blocks_per_resample = 0 !! ceiling(n_windows/block_length)
        integer :: seed = 0 !! seed of local, reproducible Park-Miller RNG
        logical :: circular = .false. !! wrap sampled blocks at the end of the window sequence
        logical :: automatic_block_length = .false. !! true if omitted length was chosen from pilot autocorrelations
        integer :: selection_lags = 0 !! number of lags retained by automatic pilot
        character(len=12) :: selection_method = 'manual' !! 'manual', 'spectral' or 'influence'
        character(len=12) :: interval_method = 'percentile' !! percentile, basic, bc, bca, studentized
        integer :: selection_mode = 0 !! influence target mode, 0 means aggregate group coherence
        integer :: selection_frequency = 0 !! influence target frequency; 0 means mean over bins
        real(dp), allocatable :: selection_influence(:) !! influence scores if influence pilot used
        integer :: status = partial_coherence_invalid_input !! success, invalid input or singular
    end type canonical_moving_block_bootstrap_result

    interface multitaper_canonical_moving_block_bootstrap
        module procedure multitaper_bootstrap_real
        module procedure multitaper_bootstrap_complex
    end interface

    public :: canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap

contains

    ! Integer arithmetic cannot overflow int64: 16807 * (2**31-2) < 2**46.
    ! The RNG is local to this call; it does not modify Fortran RANDOM_NUMBER state.
    integer function random_index(state,count) result(index)
        integer(int64), intent(inout) :: state !! mutable local Park-Miller RNG state
        integer, intent(in) :: count !! number of allowed start locations, >=1

        state = modulo(rng_multiplier*state,rng_modulus)
        index = int(modulo(state-1_int64,int(count,int64))) + 1
    end function random_index

    recursive subroutine sort_values(a,lo,hi)
        real(dp), intent(inout) :: a(:) !! finite values to sort in ascending order
        integer, intent(in) :: lo !! start index of this subrange
        integer, intent(in) :: hi !! end index of this subrange
        integer :: i, j
        real(dp) :: pivot, temp

        if (lo >= hi) return
        pivot = a(lo+(hi-lo)/2)
        i = lo
        j = hi
        do
            do while (a(i) < pivot)
                i = i+1
            end do
            do while (a(j) > pivot)
                j = j-1
            end do
            if (i > j) exit
            temp = a(i)
            a(i) = a(j)
            a(j) = temp
            i = i+1
            j = j-1
            if (i > j) exit
        end do
        if (lo < j) call sort_values(a,lo,j)
        if (i < hi) call sort_values(a,i,hi)
    end subroutine sort_values

    pure real(dp) function empirical_quantile(sorted,p) result(q)
        real(dp), intent(in) :: sorted(:) !! nonempty sorted sample
        real(dp), intent(in) :: p !! quantile probability, 0 <= p <= 1
        real(dp) :: position, weight
        integer :: lower

        position = p*real(size(sorted)-1,dp)
        lower = int(position)+1
        if (lower == size(sorted)) then
            q = sorted(lower)
        else
            weight = position-real(lower-1,dp)
            q = (1.0_dp-weight)*sorted(lower)+weight*sorted(lower+1)
        end if
    end function empirical_quantile

    ! Bisection avoids depending on a particular inverse-normal implementation.
    ! It is only used for bias-corrected percentile endpoints, never the RNG.
    pure real(dp) function normal_cdf(z) result(p)
        real(dp), intent(in) :: z !! real standard-normal argument

        p = 0.5_dp*erfc(-z/sqrt(2.0_dp))
    end function normal_cdf

    pure real(dp) function normal_quantile(p) result(z)
        real(dp), intent(in) :: p !! probability strictly between zero and one
        real(dp) :: lo, hi, mid
        integer :: step

        lo = -12.0_dp
        hi = 12.0_dp
        do step = 1, 64
            mid = 0.5_dp*(lo+hi)
            if (normal_cdf(mid) < p) then
                lo = mid
            else
                hi = mid
            end if
        end do
        z = 0.5_dp*(lo+hi)
    end function normal_quantile

    ! The caller supplies sorted bootstrap refits. Undefined BC bias offsets
    ! are reported as invalid rather than replaced by extreme endpoint bounds.
    pure subroutine interval_bounds(sorted,estimate,method,lower,upper,valid,acceleration)
        real(dp), intent(in) :: sorted(:) !! ascending finite bootstrap replicates
        real(dp), intent(in) :: estimate !! full-data estimate, expected in [0,1]
        character(*), intent(in) :: method !! percentile, basic, bc or bca
        real(dp), intent(in), optional :: acceleration !! finite block-jackknife skewness, required for bca
        real(dp), intent(out) :: lower !! pointwise lower confidence bound or NaN
        real(dp), intent(out) :: upper !! pointwise upper confidence bound or NaN
        logical, intent(out) :: valid !! whether the interval is identifiable
        real(dp) :: bias_fraction, z0, lower_probability, upper_probability, tie_tolerance
        real(dp) :: direction_low, direction_high, denominator_low, denominator_high
        integer :: nless, nequal

        lower = ieee_value(0.0_dp,ieee_quiet_nan)
        upper = lower
        valid = .true.
        select case (method)
        case ('percentile')
            lower = empirical_quantile(sorted,0.025_dp)
            upper = empirical_quantile(sorted,0.975_dp)
        case ('basic')
            lower = max(0.0_dp,2.0_dp*estimate-empirical_quantile(sorted,0.975_dp))
            upper = min(1.0_dp,2.0_dp*estimate-empirical_quantile(sorted,0.025_dp))
            ! Truncation to the parameter space may collapse an interval.
            lower = min(1.0_dp,lower)
            upper = max(0.0_dp,upper)
            if (lower > upper) then
                lower = ieee_value(0.0_dp,ieee_quiet_nan)
                upper = lower
                valid = .false.
            end if
        case ('bc','bca')
            ! Floating-point means can differ by roundoff when a resample
            ! reproduces the full spectrum. Treat those numerical ties alike.
            tie_tolerance = 64.0_dp*epsilon(estimate)*max(1.0_dp,abs(estimate))
            nless = count(sorted < estimate-tie_tolerance)
            nequal = count(abs(sorted-estimate) <= tie_tolerance)
            bias_fraction = (real(nless,dp)+0.5_dp*real(nequal,dp))/real(size(sorted),dp)
            if (bias_fraction <= 0.0_dp .or. bias_fraction >= 1.0_dp) then
                valid = .false.
                return
            end if
            z0 = normal_quantile(bias_fraction)
            if (method == 'bca') then
                if (.not. present(acceleration)) then
                    valid = .false.
                    return
                end if
                if (.not. ieee_is_finite(acceleration)) then
                    valid = .false.
                    return
                end if
                direction_low = z0-normal_975
                direction_high = z0+normal_975
                denominator_low = 1.0_dp-acceleration*direction_low
                denominator_high = 1.0_dp-acceleration*direction_high
                ! Crossing the acceleration pole invalidates the BCa map.
                if (denominator_low <= 0.0_dp .or. denominator_high <= 0.0_dp) then
                    valid = .false.
                    return
                end if
                lower_probability = normal_cdf(z0+direction_low/denominator_low)
                upper_probability = normal_cdf(z0+direction_high/denominator_high)
            else
                lower_probability = normal_cdf(2.0_dp*z0-normal_975)
                upper_probability = normal_cdf(2.0_dp*z0+normal_975)
            end if
            if (.not. ieee_is_finite(lower_probability) .or. &
                .not. ieee_is_finite(upper_probability)) then
                valid = .false.
                return
            end if
            if (method == 'bca' .and. (lower_probability <= 0.0_dp .or. &
                upper_probability >= 1.0_dp .or. lower_probability >= upper_probability)) then
                valid = .false.
                return
            end if
            lower = empirical_quantile(sorted,lower_probability)
            upper = empirical_quantile(sorted,upper_probability)
        case default
            valid = .false.
        end select
    end subroutine interval_bounds

    ! A cluster delete-one-block jackknife estimates the variance of a scalar
    ! canonical statistic. All windows in a block are removed together. For
    ! partial final blocks, weighted delete-group pseudovalues are used:
    !   z_b = [N*theta - (N-n_b)*theta_{(-b)}]/n_b,
    !   var = sum(w_b**2*(z_b-zbar)**2)/(1-sum(w_b**2)), w_b=n_b/N.
    ! The equal-block case reduces exactly to the ordinary jackknife formula.
    ! Returned NaNs mean singular/inseparable modes or unreliable variance.
    subroutine block_studentizing_errors(point,total,block_sums,block_sizes,frequencies, &
            group_a,group_b,conditioned_on,ridge,shrinkage,gap_tolerance,mode_se,group_se)
        type(canonical_coherence_result), intent(in) :: point !! estimate from all N windows
        complex(dp), intent(in) :: total(:,:,:) !! total cross-spectrum across all windows (frequency,channel,channel)
        complex(dp), intent(in) :: block_sums(:,:,:,:) !! CSD sums for each contiguous block
        integer, intent(in) :: block_sizes(:) !! number of windows in each block
        real(dp), intent(in) :: frequencies(:) !! original or band-integrated frequencies
        integer, intent(in) :: group_a(:) !! group A channel indices
        integer, intent(in) :: group_b(:) !! group B channel indices
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        real(dp), intent(in), optional :: ridge !! regularization used in point estimate
        real(dp), intent(in), optional :: shrinkage !! fixed shrinkage used in point estimate
        real(dp), intent(in) :: gap_tolerance !! isolation threshold for canonical roots
        real(dp), intent(out) :: mode_se(:,:) !! positive block-jackknife SE per (frequency,mode), otherwise NaN
        real(dp), intent(out) :: group_se(:) !! positive block-jackknife SE of aggregate, otherwise NaN
        type(multitaper_matrix_result) :: omit_spectrum
        type(canonical_coherence_result) :: leave
        real(dp), allocatable :: pseudo_modes(:,:,:), pseudo_groups(:,:), weights(:)
        real(dp) :: nan, mean_value, variance, denominator, minimum_se
        integer :: nf, nm, nc, nblocks, n, k, f, m, remaining

        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        mode_se = nan
        group_se = nan
        nf = size(point%squared_modes,1)
        nm = size(point%squared_modes,2)
        nc = size(total,2)
        nblocks = size(block_sizes)
        if (nblocks < 4 .or. size(block_sums,4) /= nblocks) return
        n = sum(block_sizes)
        if (n <= 0 .or. any(block_sizes <= 0) .or. any(block_sizes >= n)) return
        weights = real(block_sizes,dp)/real(n,dp)
        denominator = 1.0_dp-sum(weights**2)
        if (denominator <= 0.0_dp) return
        minimum_se = 64.0_dp*epsilon(1.0_dp)
        allocate(pseudo_modes(nf,nm,nblocks),pseudo_groups(nf,nblocks))
        pseudo_modes = nan
        pseudo_groups = nan
        allocate(omit_spectrum%frequencies(nf),omit_spectrum%power(nf,nc,nc))
        omit_spectrum%frequencies = frequencies
        omit_spectrum%status = partial_coherence_success
        do k = 1, nblocks
            remaining = n-block_sizes(k)
            omit_spectrum%power = (total-block_sums(:,:,:,k))/real(remaining,dp)
            leave = canonical_coherence_matrix(omit_spectrum,group_a,group_b, &
                conditioned_on,ridge=ridge,shrinkage=shrinkage)
            if (leave%status == partial_coherence_invalid_input) return
            do f = 1, nf
                if (.not. point%valid(f) .or. .not. leave%valid(f)) cycle
                if (point%coherence(f) <= gap_tolerance .or. &
                    point%coherence(f) >= 1.0_dp-gap_tolerance .or. &
                    leave%coherence(f) <= gap_tolerance .or. &
                    leave%coherence(f) >= 1.0_dp-gap_tolerance) cycle
                pseudo_groups(f,k) = (real(n,dp)*point%coherence(f)- &
                    real(remaining,dp)*leave%coherence(f))/real(block_sizes(k),dp)
                do m = 1, nm
                    if (.not. isolated_mode(point%squared_modes(f,:),m,gap_tolerance) .or. &
                        .not. isolated_mode(leave%squared_modes(f,:),m,gap_tolerance)) cycle
                    pseudo_modes(f,m,k) = (real(n,dp)*point%squared_modes(f,m)- &
                        real(remaining,dp)*leave%squared_modes(f,m))/real(block_sizes(k),dp)
                end do
            end do
        end do
        do f = 1, nf
            if (all(ieee_is_finite(pseudo_groups(f,:)))) then
                mean_value = sum(weights*pseudo_groups(f,:))
                variance = sum(weights**2*(pseudo_groups(f,:)-mean_value)**2)/denominator
                if (ieee_is_finite(variance) .and. variance > minimum_se**2) group_se(f) = sqrt(variance)
            end if
            do m = 1, nm
                if (.not. all(ieee_is_finite(pseudo_modes(f,m,:)))) cycle
                mean_value = sum(weights*pseudo_modes(f,m,:))
                variance = sum(weights**2*(pseudo_modes(f,m,:)-mean_value)**2)/denominator
                if (ieee_is_finite(variance) .and. variance > minimum_se**2) mode_se(f,m) = sqrt(variance)
            end do
        end do
    end subroutine block_studentizing_errors

    ! Weighted, centered delete-group pseudovalue contributions. With
    ! equal-size groups this is the usual jackknife BCa acceleration.
    ! Unequal final blocks give a first-order cluster-influence approximation.
    pure real(dp) function block_acceleration(pseudo,weights) result(acceleration)
        real(dp), intent(in) :: pseudo(:) !! finite delete-group pseudovalues
        real(dp), intent(in) :: weights(:) !! corresponding positive fractions summing to one
        real(dp) :: contributions(size(pseudo)), scale, second, centered

        acceleration = ieee_value(0.0_dp,ieee_quiet_nan)
        if (.not. all(ieee_is_finite(pseudo))) return
        centered = sum(weights*pseudo)
        contributions = weights*(pseudo-centered)
        scale = maxval(abs(contributions))
        if (.not. ieee_is_finite(scale) .or. scale <= 64.0_dp*epsilon(1.0_dp)) return
        contributions = contributions/scale
        second = sum(contributions**2)
        if (second <= 0.0_dp) return
        acceleration = sum(contributions**3)/(6.0_dp*second**1.5_dp)
    end function block_acceleration

    ! Acceleration is estimated once from the ORIGINAL complete windows.
    ! It is not re-estimated for each moving-block bootstrap replicate.
    ! Inseparable roots, boundary coherences and failed deletion refits
    ! have NaN accelerations and therefore no BCa interval.
    subroutine block_bca_accelerations(point,total,block_sums,block_sizes,frequencies, &
            group_a,group_b,conditioned_on,ridge,shrinkage,gap_tolerance,mode_accel,group_accel)
        type(canonical_coherence_result), intent(in) :: point !! full-window estimator
        complex(dp), intent(in) :: total(:,:,:) !! sum of complete-window CSD matrices
        complex(dp), intent(in) :: block_sums(:,:,:,:) !! original consecutive block sums
        integer, intent(in) :: block_sizes(:) !! original block lengths, including final remainder
        real(dp), intent(in) :: frequencies(:) !! frequencies or integrated band centers
        integer, intent(in) :: group_a(:) !! first channel group indices
        integer, intent(in) :: group_b(:) !! second channel group indices
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels
        real(dp), intent(in), optional :: ridge !! fixed ridge applied in all refits
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        real(dp), intent(in) :: gap_tolerance !! minimum interior and eigenvalue gap
        real(dp), intent(out) :: mode_accel(:,:) !! acceleration for each squared mode or NaN
        real(dp), intent(out) :: group_accel(:) !! acceleration for aggregate coherence or NaN
        type(multitaper_matrix_result) :: omit_spectrum
        type(canonical_coherence_result) :: leave
        real(dp), allocatable :: pseudo_modes(:,:,:), pseudo_groups(:,:), weights(:)
        real(dp) :: nan
        integer :: nf, nm, nc, nblocks, n, k, f, m, remaining

        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        mode_accel = nan
        group_accel = nan
        nf = size(point%squared_modes,1)
        nm = size(point%squared_modes,2)
        nc = size(total,2)
        nblocks = size(block_sizes)
        if (nblocks < 4 .or. size(block_sums,4) /= nblocks) return
        n = sum(block_sizes)
        if (n <= 0 .or. any(block_sizes <= 0) .or. any(block_sizes >= n)) return
        weights = real(block_sizes,dp)/real(n,dp)
        allocate(pseudo_modes(nf,nm,nblocks),pseudo_groups(nf,nblocks))
        pseudo_modes = nan
        pseudo_groups = nan
        allocate(omit_spectrum%frequencies(nf),omit_spectrum%power(nf,nc,nc))
        omit_spectrum%frequencies = frequencies
        omit_spectrum%status = partial_coherence_success
        do k = 1, nblocks
            remaining = n-block_sizes(k)
            omit_spectrum%power = (total-block_sums(:,:,:,k))/real(remaining,dp)
            leave = canonical_coherence_matrix(omit_spectrum,group_a,group_b, &
                conditioned_on,ridge=ridge,shrinkage=shrinkage)
            if (leave%status == partial_coherence_invalid_input) return
            do f = 1, nf
                if (.not. point%valid(f) .or. .not. leave%valid(f)) cycle
                if (point%coherence(f) <= gap_tolerance .or. &
                    point%coherence(f) >= 1.0_dp-gap_tolerance .or. &
                    leave%coherence(f) <= gap_tolerance .or. &
                    leave%coherence(f) >= 1.0_dp-gap_tolerance) cycle
                pseudo_groups(f,k) = (real(n,dp)*point%coherence(f)- &
                    real(remaining,dp)*leave%coherence(f))/real(block_sizes(k),dp)
                do m = 1, nm
                    if (.not. isolated_mode(point%squared_modes(f,:),m,gap_tolerance) .or. &
                        .not. isolated_mode(leave%squared_modes(f,:),m,gap_tolerance)) cycle
                    pseudo_modes(f,m,k) = (real(n,dp)*point%squared_modes(f,m)- &
                        real(remaining,dp)*leave%squared_modes(f,m))/real(block_sizes(k),dp)
                end do
            end do
        end do
        do f = 1, nf
            group_accel(f) = block_acceleration(pseudo_groups(f,:),weights)
            do m = 1, nm
                mode_accel(f,m) = block_acceleration(pseudo_modes(f,m,:),weights)
            end do
        end do
    end subroutine block_bca_accelerations

    ! A bootstrap-t interval reverses the quantiles of
    ! (theta_star-theta_original)/SE_star, scaled by the original SE.
    subroutine studentized_bounds(tvalues,estimate,standard_error,lower,upper,valid)
        real(dp), intent(in) :: tvalues(:) !! one finite studentized statistic per replicate
        real(dp), intent(in) :: estimate !! full-sample canonical statistic in [0,1]
        real(dp), intent(in) :: standard_error !! original delete-one-block SE
        real(dp), intent(out) :: lower !! lower endpoint or NaN
        real(dp), intent(out) :: upper !! upper endpoint or NaN
        logical, intent(out) :: valid !! true if all studentizing SEs identifiable
        real(dp) :: nan
        real(dp), allocatable :: ordered(:)

        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        lower = nan
        upper = nan
        valid = .false.
        if (.not. ieee_is_finite(standard_error) .or. standard_error <= 0.0_dp) return
        if (.not. all(ieee_is_finite(tvalues))) return
        ordered = tvalues
        call sort_values(ordered,1,size(ordered))
        lower = max(0.0_dp,estimate-standard_error*empirical_quantile(ordered,0.975_dp))
        upper = min(1.0_dp,estimate-standard_error*empirical_quantile(ordered,0.025_dp))
        lower = min(1.0_dp,lower)
        upper = max(0.0_dp,upper)
        if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper) .or. lower > upper) then
            lower = nan
            upper = nan
            return
        end if
        valid = .true.
    end subroutine studentized_bounds

    pure logical function isolated_mode(modes,k,tolerance) result(ok)
        real(dp), intent(in) :: modes(:) !! canonical roots in descending order
        integer, intent(in) :: k !! checked mode index, 1-based
        real(dp), intent(in) :: tolerance !! minimum boundary and adjacent-root gap

        ok = .false.
        if (.not. all(ieee_is_finite(modes))) return
        if (modes(k) <= tolerance .or. modes(k) >= 1.0_dp-tolerance) return
        if (k > 1) then
            if (modes(k-1)-modes(k) <= tolerance) return
        end if
        if (k < size(modes)) then
            if (modes(k)-modes(k+1) <= tolerance) return
        end if
        ok = .true.
    end function isolated_mode

    function canonical_moving_block_bootstrap(spectra,group_a,group_b,block_length, &
            n_resamples,conditioned_on,ridge,shrinkage,bands,average,gap_tolerance,seed,circular, &
            selection_method,selection_mode,selection_frequency,interval_method) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectra !! (freq,window,channel,channel) CSD
        integer, intent(in) :: group_a(:) !! first nonempty group of channel indices, one-based
        integer, intent(in) :: group_b(:) !! second nonempty, disjoint group indices
        integer, intent(in), optional :: block_length !! optional explicit block length; omitted selects from pilot
        integer, intent(in), optional :: n_resamples !! bootstrap replicates, >=20; default 399
        integer, intent(in), optional :: conditioned_on(:) !! optional disjoint conditioning channels
        real(dp), intent(in), optional :: ridge !! fixed trace-relative ridge, >=0
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage, [0,1]
        real(dp), intent(in), optional :: bands(:,:) !! optional integration bounds, (2,nband)
        logical, intent(in), optional :: average !! return band-averaged rather than integrated CSDs
        real(dp), intent(in), optional :: gap_tolerance !! minimum isolated-root separation, default 1e-6
        integer, intent(in), optional :: seed !! deterministic independent local seed, [1,2147483646]
        logical, intent(in), optional :: circular !! allow a block to wrap around the final window (default false)
        character(*), intent(in), optional :: selection_method !! automatic selector: spectral (default) or influence
        integer, intent(in), optional :: selection_mode !! influence: 0 for group coherence, 1..nmode for root
        integer, intent(in), optional :: selection_frequency !! influence: 0 for frequency mean, else one bin
        character(*), intent(in), optional :: interval_method !! percentile (default), basic, bc, bca, studentized
        type(canonical_moving_block_bootstrap_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: work
        type(multitaper_matrix_result) :: mean_spectrum
        type(canonical_coherence_result) :: point, fit
        type(canonical_block_length_result) :: choice
        real(dp), allocatable :: roots(:,:,:), groups(:,:), sorted(:)
        real(dp), allocatable :: student_t_modes(:,:,:), student_t_groups(:,:)
        real(dp), allocatable :: replicate_mode_se(:,:), replicate_group_se(:)
        complex(dp), allocatable :: sampled_block_sums(:,:,:,:), original_block_sums(:,:,:,:)
        complex(dp), allocatable :: sampled_total(:,:,:)
        integer, allocatable :: sample_block_sizes(:), original_block_sizes(:)
        complex(dp), allocatable :: totals(:,:,:)
        real(dp) :: nan, tol, avg
        integer(int64) :: state
        integer :: nf, nw, nc, nm, nr, nstarts, nsampled, nblocks
        integer :: r, block, offset, start, window, filled, f, mode, seed_value, length
        integer :: n_orig_blocks, student_offset
        logical :: wrap, regular, automatic, interval_ok, studentized, bca
        character(len=12) :: method
        character(len=12) :: bounds_method

        bounds_method = 'percentile'
        if (present(interval_method)) then
            if (len_trim(interval_method) > len(bounds_method)) return
            bounds_method = trim(interval_method)
        end if
        if (bounds_method /= 'percentile' .and. bounds_method /= 'basic' .and. &
            bounds_method /= 'bc' .and. bounds_method /= 'bca' .and. &
            bounds_method /= 'studentized') return
        studentized = bounds_method == 'studentized'
        bca = bounds_method == 'bca'
        if (spectra%status /= partial_coherence_success) return
        if (.not. allocated(spectra%power) .or. .not. allocated(spectra%times) .or. &
            .not. allocated(spectra%frequencies)) return
        nf = size(spectra%frequencies)
        nw = size(spectra%times)
        nc = size(spectra%power,3)
        if (nf < 1 .or. nw < 4 .or. nc < 2) return
        if (size(spectra%power,1) /= nf .or. size(spectra%power,2) /= nw .or. &
            size(spectra%power,4) /= nc) return
        automatic = .not. present(block_length)
        method = 'spectral'
        if (present(selection_method)) then
            if (len_trim(selection_method) > len(method)) return
            method = trim(selection_method)
        end if
        if (method /= 'spectral' .and. method /= 'influence') return
        if (.not. automatic .and. method /= 'spectral') return
        if (method /= 'influence') then
            if (present(selection_mode) .or. present(selection_frequency)) return
        end if
        if (present(block_length)) then
            length = block_length
            if (length < 1 .or. length >= nw) return
        end if
        nr = default_resamples
        if (present(n_resamples)) nr = n_resamples
        if (nr < 20) return
        seed_value = default_seed
        if (present(seed)) seed_value = seed
        if (seed_value < 1 .or. int(seed_value,int64) >= rng_modulus) return
        tol = 1.0e-6_dp
        if (present(gap_tolerance)) tol = gap_tolerance
        if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp .or. tol >= 0.5_dp) return
        if (.not. all(ieee_is_finite(spectra%times)) .or. &
            .not. all(ieee_is_finite(spectra%frequencies))) return
        if (any(spectra%times(2:) <= spectra%times(:nw-1))) return
        if (.not. all(ieee_is_finite(real(spectra%power,dp))) .or. &
            .not. all(ieee_is_finite(aimag(spectra%power)))) return

        wrap = .false.
        if (present(circular)) wrap = circular
        if (present(bands)) then
            work = integrate_tf_spectral_bands(spectra,bands,average)
        else
            work = spectra
        end if
        if (work%status /= partial_coherence_success) return
        if (automatic) then
            if (method == 'influence') then
                choice = select_canonical_influence_block_length(work,group_a,group_b, &
                    conditioned_on=conditioned_on,ridge=ridge,shrinkage=shrinkage, &
                    mode_index=selection_mode,frequency_index=selection_frequency)
            else
                choice = select_canonical_block_length(work)
            end if
            if (choice%status /= partial_coherence_success) return
            length = choice%block_length
        end if
        nstarts = nw-length+1
        if (wrap) nstarts = nw
        nblocks = (nw-1)/length+1
        nf = size(work%frequencies)
        allocate(mean_spectrum%frequencies(nf),mean_spectrum%power(nf,nc,nc),totals(nf,nc,nc))
        mean_spectrum%frequencies = work%frequencies
        totals = sum(work%power,dim=2)
        mean_spectrum%power = totals/real(nw,dp)
        mean_spectrum%status = partial_coherence_success
        point = canonical_coherence_matrix(mean_spectrum,group_a,group_b, &
            conditioned_on,ridge=ridge,shrinkage=shrinkage)
        if (point%status == partial_coherence_invalid_input) return
        nm = size(point%squared_modes,2)
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        allocate(out%frequencies(nf),out%squared_modes(nf,nm),out%stderr(nf,nm), &
            out%bias(nf,nm),out%lower95(nf,nm),out%upper95(nf,nm),out%valid(nf,nm), &
            out%coherence(nf),out%coherence_stderr(nf),out%coherence_bias(nf), &
            out%coherence_lower95(nf),out%coherence_upper95(nf),out%coherence_valid(nf), &
            roots(nf,nm,nr),groups(nf,nr),sorted(nr))
        out%frequencies = point%frequencies
        out%squared_modes = point%squared_modes
        out%coherence = point%coherence
        out%stderr = nan
        out%bias = nan
        out%lower95 = nan
        out%upper95 = nan
        out%valid = .false.
        out%coherence_stderr = nan
        out%coherence_bias = nan
        out%coherence_lower95 = nan
        out%coherence_upper95 = nan
        out%coherence_valid = .false.
        out%n_windows = nw
        out%n_resamples = nr
        out%block_length = length
        out%automatic_block_length = automatic
        if (automatic) then
            out%selection_lags = choice%lags_used
            out%selection_method = method
            if (method == 'influence') then
                out%selection_mode = choice%mode_index
                out%selection_frequency = choice%frequency_index
                out%selection_influence = choice%influence
            end if
        end if
        out%blocks_per_resample = nblocks
        out%seed = seed_value
        out%interval_method = bounds_method
        out%circular = wrap
        out%status = point%status
        roots = nan
        groups = nan

        if (studentized .or. bca) then
            if (studentized) then
                allocate(out%studentizing_stderr(nf,nm),out%coherence_studentizing_stderr(nf), &
                    student_t_modes(nf,nm,nr),student_t_groups(nf,nr), &
                    replicate_mode_se(nf,nm),replicate_group_se(nf), &
                    sampled_block_sums(nf,nc,nc,nblocks),sample_block_sizes(nblocks), &
                    sampled_total(nf,nc,nc))
                student_t_modes = nan
                student_t_groups = nan
            else
                allocate(out%bca_acceleration(nf,nm),out%coherence_bca_acceleration(nf))
            end if
            ! Original data are partitioned into contiguous complete-window
            ! blocks, with a possibly shorter final block.
            n_orig_blocks = nblocks
            allocate(original_block_sums(nf,nc,nc,n_orig_blocks), &
                original_block_sizes(n_orig_blocks))
            student_offset = 0
            do block = 1, n_orig_blocks
                original_block_sizes(block) = min(length,nw-student_offset)
                original_block_sums(:,:,:,block) = &
                    sum(work%power(:,student_offset+1:student_offset+original_block_sizes(block),:,:),dim=2)
                student_offset = student_offset+original_block_sizes(block)
            end do
            if (studentized) then
                call block_studentizing_errors(point,totals,original_block_sums,original_block_sizes, &
                    work%frequencies,group_a,group_b,conditioned_on,ridge,shrinkage,tol, &
                    out%studentizing_stderr,out%coherence_studentizing_stderr)
            else
                call block_bca_accelerations(point,totals,original_block_sums,original_block_sizes, &
                    work%frequencies,group_a,group_b,conditioned_on,ridge,shrinkage,tol, &
                    out%bca_acceleration,out%coherence_bca_acceleration)
            end if
        end if

        state = int(seed_value,int64)
        do r = 1, nr
            mean_spectrum%power = (0.0_dp,0.0_dp)
            filled = 0
            do block = 1, nblocks
                start = random_index(state,nstarts)
                nsampled = min(length,nw-filled)
                if (studentized) then
                    sample_block_sizes(block) = nsampled
                    sampled_block_sums(:,:,:,block) = (0.0_dp,0.0_dp)
                end if
                do offset = 0, nsampled-1
                    window = start+offset
                    if (wrap) window = modulo(window-1,nw)+1
                    mean_spectrum%power = mean_spectrum%power+work%power(:,window,:,:)
                    if (studentized) sampled_block_sums(:,:,:,block) = &
                        sampled_block_sums(:,:,:,block)+work%power(:,window,:,:)
                end do
                filled = filled+nsampled
            end do
            if (studentized) sampled_total = mean_spectrum%power
            mean_spectrum%power = mean_spectrum%power/real(nw,dp)
            fit = canonical_coherence_matrix(mean_spectrum,group_a,group_b, &
                conditioned_on,ridge=ridge,shrinkage=shrinkage)
            if (fit%status == partial_coherence_invalid_input) then
                out%status = partial_coherence_invalid_input
                return
            end if
            if (studentized) then
                call block_studentizing_errors(fit,sampled_total,sampled_block_sums,sample_block_sizes, &
                    work%frequencies,group_a,group_b,conditioned_on,ridge,shrinkage,tol, &
                    replicate_mode_se,replicate_group_se)
            end if
            do f = 1, nf
                if (.not. fit%valid(f)) cycle
                roots(f,:,r) = fit%squared_modes(f,:)
                groups(f,r) = fit%coherence(f)
                if (.not. studentized) cycle
                if (ieee_is_finite(replicate_group_se(f))) student_t_groups(f,r) = &
                    (fit%coherence(f)-point%coherence(f))/replicate_group_se(f)
                do mode = 1, nm
                    if (.not. ieee_is_finite(replicate_mode_se(f,mode))) cycle
                    student_t_modes(f,mode,r) = &
                        (fit%squared_modes(f,mode)-point%squared_modes(f,mode))/replicate_mode_se(f,mode)
                end do
            end do
        end do

        do f = 1, nf
            if (.not. point%valid(f)) cycle
            if (.not. all(ieee_is_finite(roots(f,:,:)))) cycle
            if (.not. all(ieee_is_finite(groups(f,:)))) cycle
            sorted = groups(f,:)
            call sort_values(sorted,1,nr)
            avg = sum(groups(f,:))/real(nr,dp)
            out%coherence_bias(f) = avg-point%coherence(f)
            out%coherence_stderr(f) = sqrt(sum((groups(f,:)-avg)**2)/real(nr-1,dp))
            if (studentized) then
                call studentized_bounds(student_t_groups(f,:),point%coherence(f), &
                    out%coherence_studentizing_stderr(f),out%coherence_lower95(f), &
                    out%coherence_upper95(f),interval_ok)
            else if (bca) then
                call interval_bounds(sorted,point%coherence(f),bounds_method, &
                    out%coherence_lower95(f),out%coherence_upper95(f),interval_ok, &
                    out%coherence_bca_acceleration(f))
            else
                call interval_bounds(sorted,point%coherence(f),bounds_method, &
                    out%coherence_lower95(f),out%coherence_upper95(f),interval_ok)
            end if
            out%coherence_valid(f) = interval_ok
            do mode = 1, nm
                regular = isolated_mode(point%squared_modes(f,:),mode,tol)
                if (.not. regular) cycle
                do r = 1, nr
                    if (.not. isolated_mode(roots(f,:,r),mode,tol)) then
                        regular = .false.
                        exit
                    end if
                end do
                if (.not. regular) cycle
                sorted = roots(f,mode,:)
                call sort_values(sorted,1,nr)
                avg = sum(roots(f,mode,:))/real(nr,dp)
                out%bias(f,mode) = avg-point%squared_modes(f,mode)
                out%stderr(f,mode) = sqrt(sum((roots(f,mode,:)-avg)**2)/real(nr-1,dp))
                if (studentized) then
                    call studentized_bounds(student_t_modes(f,mode,:),point%squared_modes(f,mode), &
                        out%studentizing_stderr(f,mode),out%lower95(f,mode), &
                        out%upper95(f,mode),interval_ok)
                else if (bca) then
                    call interval_bounds(sorted,point%squared_modes(f,mode),bounds_method, &
                        out%lower95(f,mode),out%upper95(f,mode),interval_ok, &
                        out%bca_acceleration(f,mode))
                else
                    call interval_bounds(sorted,point%squared_modes(f,mode),bounds_method, &
                        out%lower95(f,mode),out%upper95(f,mode),interval_ok)
                end if
                out%valid(f,mode) = interval_ok
            end do
        end do
        if (.not. all(out%coherence_valid)) out%status = partial_coherence_singular
    end function canonical_moving_block_bootstrap

    function multitaper_bootstrap_real(x,nw,nperseg,group_a,group_b,block_length,n_resamples, &
            conditioned_on,kmax,fs,nfft,detrend,weighting,axis,detrend_mode, &
            ridge,shrinkage,bands,average,gap_tolerance,seed,circular, &
            selection_method,selection_mode,selection_frequency,interval_method) result(out)
        real(dp), intent(in) :: x(:,:) !! input real samples, axis-selected samples and channels
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples per nonoverlapping spectral-estimation window
        integer, intent(in) :: group_a(:) !! 1-based group A channel indices
        integer, intent(in) :: group_b(:) !! 1-based group B channel indices
        integer, intent(in), optional :: block_length !! optional explicit length; omitted selects automatically
        integer, intent(in), optional :: n_resamples !! number of bootstrap replicates, default 399
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels
        integer, intent(in), optional :: kmax !! DPSS taper count
        real(dp), intent(in), optional :: fs !! sample frequency, Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove window mean
        character(*), intent(in), optional :: weighting !! equal, eigen or adaptive
        integer, intent(in), optional :: axis !! sample axis (1 or 2)
        character(*), intent(in), optional :: detrend_mode !! none, constant or linear
        real(dp), intent(in), optional :: ridge !! fixed trace-relative ridge
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        real(dp), intent(in), optional :: bands(:,:) !! optional frequency bands, shape (2,nband)
        logical, intent(in), optional :: average !! average band CSD instead of integral
        real(dp), intent(in), optional :: gap_tolerance !! root isolation threshold
        integer, intent(in), optional :: seed !! local reproducible RNG seed
        logical, intent(in), optional :: circular !! allow blocks to wrap around sequence
        character(*), intent(in), optional :: selection_method !! automatic spectral/influence pilot
        integer, intent(in), optional :: selection_mode !! 0 for group coherence; positive for a mode
        integer, intent(in), optional :: selection_frequency !! 0 for frequency mean; positive for one bin
        character(*), intent(in), optional :: interval_method !! percentile (default), basic, bc, bca, studentized
        type(canonical_moving_block_bootstrap_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectra

        spectra = multitaper_spectral_matrix_spectrogram(x,nw,nperseg, &
            noverlap=0,kmax=kmax,fs=fs,nfft=nfft,detrend=detrend,weighting=weighting, &
            axis=axis,detrend_mode=detrend_mode)
        out = canonical_moving_block_bootstrap(spectra,group_a,group_b,block_length,n_resamples, &
            conditioned_on,ridge,shrinkage,bands,average,gap_tolerance,seed,circular, &
            selection_method,selection_mode,selection_frequency,interval_method)
    end function multitaper_bootstrap_real

    function multitaper_bootstrap_complex(x,nw,nperseg,group_a,group_b,block_length,n_resamples, &
            conditioned_on,kmax,fs,nfft,detrend,weighting,axis,detrend_mode, &
            ridge,shrinkage,bands,average,gap_tolerance,seed,circular, &
            selection_method,selection_mode,selection_frequency,interval_method) result(out)
        complex(dp), intent(in) :: x(:,:) !! input complex samples, axis-selected samples and channels
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples per nonoverlapping spectral-estimation window
        integer, intent(in) :: group_a(:) !! 1-based group A channel indices
        integer, intent(in) :: group_b(:) !! 1-based group B channel indices
        integer, intent(in), optional :: block_length !! optional explicit length; omitted selects automatically
        integer, intent(in), optional :: n_resamples !! number of bootstrap replicates, default 399
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels
        integer, intent(in), optional :: kmax !! DPSS taper count
        real(dp), intent(in), optional :: fs !! sample frequency, Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove window mean
        character(*), intent(in), optional :: weighting !! equal, eigen or adaptive
        integer, intent(in), optional :: axis !! sample axis (1 or 2)
        character(*), intent(in), optional :: detrend_mode !! none, constant or linear
        real(dp), intent(in), optional :: ridge !! fixed trace-relative ridge
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        real(dp), intent(in), optional :: bands(:,:) !! optional frequency bands, shape (2,nband)
        logical, intent(in), optional :: average !! average band CSD instead of integral
        real(dp), intent(in), optional :: gap_tolerance !! root isolation threshold
        integer, intent(in), optional :: seed !! local reproducible RNG seed
        logical, intent(in), optional :: circular !! allow blocks to wrap around sequence
        character(*), intent(in), optional :: selection_method !! automatic spectral/influence pilot
        integer, intent(in), optional :: selection_mode !! 0 for group coherence; positive for a mode
        integer, intent(in), optional :: selection_frequency !! 0 for frequency mean; positive for one bin
        character(*), intent(in), optional :: interval_method !! percentile (default), basic, bc, bca, studentized
        type(canonical_moving_block_bootstrap_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectra

        spectra = multitaper_spectral_matrix_spectrogram(x,nw,nperseg, &
            noverlap=0,kmax=kmax,fs=fs,nfft=nfft,detrend=detrend,weighting=weighting, &
            axis=axis,detrend_mode=detrend_mode)
        out = canonical_moving_block_bootstrap(spectra,group_a,group_b,block_length,n_resamples, &
            conditioned_on,ridge,shrinkage,bands,average,gap_tolerance,seed,circular, &
            selection_method,selection_mode,selection_frequency,interval_method)
    end function multitaper_bootstrap_complex

end module scifort_signal_canonical_bootstrap
