! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Exact stationary VAR(p) population targets for a fixed-weight multitaper
! spectral matrix. The expectation is taken BEFORE canonical analysis.
! This is a SciFort extension, not a scipy.signal public API.
!
! Weights are fixed independently of the data and frequency. Consequently
! the band integrals of the expected spectral matrix remain analytic.
! Adaptive/data-dependent taper weighting is deliberately not represented.
module scifort_signal_varp_multitaper_target
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_matrix_result
    use scifort_signal_varp_taper_spectrum, only : varp_taper_spectrum_result, &
        varp_expected_taper_spectrum, varp_expected_band_taper_spectrum, &
        varp_spectrum_success
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, &
        canonical_coherence_matrix
    implicit none
    private

    public :: varp_expected_multitaper_spectrum
    public :: varp_expected_band_multitaper_spectrum
    public :: varp_expected_multitaper_canonical_coherence
    public :: varp_expected_band_multitaper_canonical_coherence

contains

    ! Dividing by the largest supplied weight prevents overflow of the weight
    ! sum, while preserving relative weights including exact zero weights.
    pure subroutine checked_weights(tapers,taper_weights,weights,ok)
        real(dp), intent(in) :: tapers(:,:) !! samples by taper, at least one of each
        real(dp), intent(in), optional :: taper_weights(:) !! fixed nonnegative taper weights
        real(dp), allocatable, intent(out) :: weights(:) !! positive weights summing to one
        logical, intent(out) :: ok !! true when bank and weights are valid
        real(dp) :: scale, total
        integer :: k

        ok = .false.
        if (size(tapers,1) < 1 .or. size(tapers,2) < 1) return
        if (.not. all(ieee_is_finite(tapers))) return
        allocate(weights(size(tapers,2)))
        weights = 1.0_dp
        if (present(taper_weights)) then
            if (size(taper_weights) /= size(weights)) return
            if (.not. all(ieee_is_finite(taper_weights))) return
            if (any(taper_weights < 0.0_dp)) return
            weights = taper_weights
        end if
        scale = maxval(weights)
        if (scale <= 0.0_dp) return
        weights = weights/scale
        total = sum(weights)
        if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) return
        weights = weights/total
        ! Zero-weight windows are excluded from the estimator. Every selected
        ! taper must have positive, finite energy for either normalization.
        do k = 1, size(weights)
            if (weights(k) == 0.0_dp) cycle
            if (maxval(abs(tapers(:,k))) == 0.0_dp) return
            if (.not. ieee_is_finite(sum(tapers(:,k)**2))) return
            if (sum(tapers(:,k)**2) <= 0.0_dp) return
        end do
        ok = .true.
    end subroutine checked_weights

    function varp_expected_multitaper_spectrum(a,innovation_covariance,tapers,frequencies, &
            taper_weights,normalized) result(out)
        real(dp), intent(in) :: a(:,:,:) !! stable VAR(p) coefficient matrices (n,n,p)
        real(dp), intent(in) :: innovation_covariance(:,:) !! positive-definite innovation matrix
        real(dp), intent(in) :: tapers(:,:) !! deterministic bank (window samples, taper number)
        real(dp), intent(in) :: frequencies(:) !! output frequencies, in cycles/sample
        real(dp), intent(in), optional :: taper_weights(:) !! nonnegative fixed weights, default equal
        logical, intent(in), optional :: normalized !! normalize energy of EACH taper, default false
        type(varp_taper_spectrum_result) :: out
        type(varp_taper_spectrum_result) :: single, working
        real(dp), allocatable :: weights(:)
        integer :: k
        logical :: ok

        call checked_weights(tapers,taper_weights,weights,ok)
        if (.not. ok) return
        ! Computing each finite-window target with the existing exact solver
        ! also retains its original stationarity and covariance validation.
        do k = 1, size(weights)
            if (weights(k) == 0.0_dp) cycle
            single = varp_expected_taper_spectrum(a,innovation_covariance,tapers(:,k), &
                frequencies,normalized)
            if (single%status /= varp_spectrum_success) then
                out%status = single%status
                return
            end if
            if (.not. allocated(working%power)) then
                working%frequencies = single%frequencies
                working%lag0_covariance = single%lag0_covariance
                allocate(working%power, mold=single%power)
                working%power = (0.0_dp,0.0_dp)
            end if
            working%power = working%power + weights(k)*single%power
        end do
        working%status = varp_spectrum_success
        out = working
    end function varp_expected_multitaper_spectrum

    function varp_expected_band_multitaper_spectrum(a,innovation_covariance,tapers,bands, &
            taper_weights,average,normalized) result(out)
        real(dp), intent(in) :: a(:,:,:) !! stable VAR(p) coefficient matrices (n,n,p)
        real(dp), intent(in) :: innovation_covariance(:,:) !! positive-definite innovation matrix
        real(dp), intent(in) :: tapers(:,:) !! deterministic bank (window samples, taper number)
        real(dp), intent(in) :: bands(:,:) !! band lower/upper edges (2,nband), cycles/sample
        real(dp), intent(in), optional :: taper_weights(:) !! fixed nonnegative weights, default equal
        logical, intent(in), optional :: average !! divide band integral by width; default false
        logical, intent(in), optional :: normalized !! divide EACH taper by its energy; default false
        type(varp_taper_spectrum_result) :: out
        type(varp_taper_spectrum_result) :: single, working
        real(dp), allocatable :: weights(:)
        integer :: k
        logical :: ok

        call checked_weights(tapers,taper_weights,weights,ok)
        if (.not. ok) return
        do k = 1, size(weights)
            if (weights(k) == 0.0_dp) cycle
            single = varp_expected_band_taper_spectrum(a,innovation_covariance,tapers(:,k), &
                bands,average=average,normalized=normalized)
            if (single%status /= varp_spectrum_success) then
                out%status = single%status
                return
            end if
            if (.not. allocated(working%power)) then
                working%frequencies = single%frequencies
                working%bands = single%bands
                working%lag0_covariance = single%lag0_covariance
                allocate(working%power, mold=single%power)
                working%power = (0.0_dp,0.0_dp)
            end if
            working%power = working%power + weights(k)*single%power
        end do
        working%status = varp_spectrum_success
        out = working
    end function varp_expected_band_multitaper_spectrum

    function varp_expected_multitaper_canonical_coherence(a,innovation_covariance, &
            tapers,frequencies,group_a,group_b,conditioned_on,taper_weights,ridge, &
            shrinkage,normalized,return_vectors) result(out)
        real(dp), intent(in) :: a(:,:,:) !! stable VAR(p) coefficients (n,n,p)
        real(dp), intent(in) :: innovation_covariance(:,:) !! positive-definite innovation matrix
        real(dp), intent(in) :: tapers(:,:) !! fixed taper bank (window samples, taper)
        real(dp), intent(in) :: frequencies(:) !! output frequencies, cycles/sample
        integer, intent(in) :: group_a(:) !! one-based indices of first channel group
        integer, intent(in) :: group_b(:) !! one-based disjoint second channel group
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels
        real(dp), intent(in), optional :: taper_weights(:) !! fixed nonnegative weights
        real(dp), intent(in), optional :: ridge !! nonnegative relative ridge loading
        real(dp), intent(in), optional :: shrinkage !! diagonal shrinkage in [0,1]
        logical, intent(in), optional :: normalized !! energy-normalize each taper; default false
        logical, intent(in), optional :: return_vectors !! return canonical directions; default false
        type(canonical_coherence_result) :: out
        type(varp_taper_spectrum_result) :: exact
        type(multitaper_matrix_result) :: spectrum

        exact = varp_expected_multitaper_spectrum(a,innovation_covariance,tapers, &
            frequencies,taper_weights=taper_weights,normalized=normalized)
        if (exact%status /= varp_spectrum_success) return
        spectrum%frequencies = exact%frequencies
        spectrum%power = exact%power
        spectrum%status = spectral_status_success
        out = canonical_coherence_matrix(spectrum,group_a,group_b, &
            conditioned_on=conditioned_on,ridge=ridge,shrinkage=shrinkage, &
            return_vectors=return_vectors)
    end function varp_expected_multitaper_canonical_coherence

    function varp_expected_band_multitaper_canonical_coherence(a,innovation_covariance, &
            tapers,bands,group_a,group_b,conditioned_on,taper_weights,average, &
            ridge,shrinkage,normalized,return_vectors) result(out)
        real(dp), intent(in) :: a(:,:,:) !! stable VAR(p) coefficients (n,n,p)
        real(dp), intent(in) :: innovation_covariance(:,:) !! positive-definite innovation matrix
        real(dp), intent(in) :: tapers(:,:) !! fixed taper bank (window samples, taper)
        real(dp), intent(in) :: bands(:,:) !! exact boundaries (2,nband), cycles/sample
        integer, intent(in) :: group_a(:) !! one-based indices of first channel group
        integer, intent(in) :: group_b(:) !! one-based disjoint second channel group
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels
        real(dp), intent(in), optional :: taper_weights(:) !! fixed nonnegative weights
        logical, intent(in), optional :: average !! return mean rather than band integral
        real(dp), intent(in), optional :: ridge !! nonnegative relative ridge loading
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage in [0,1]
        logical, intent(in), optional :: normalized !! energy-normalize each taper; default false
        logical, intent(in), optional :: return_vectors !! return canonical directions; default false
        type(canonical_coherence_result) :: out
        type(varp_taper_spectrum_result) :: exact
        type(multitaper_matrix_result) :: spectrum

        exact = varp_expected_band_multitaper_spectrum(a,innovation_covariance,tapers,bands, &
            taper_weights=taper_weights,average=average,normalized=normalized)
        if (exact%status /= varp_spectrum_success) return
        spectrum%frequencies = exact%frequencies
        spectrum%power = exact%power
        spectrum%status = spectral_status_success
        out = canonical_coherence_matrix(spectrum,group_a,group_b, &
            conditioned_on=conditioned_on,ridge=ridge,shrinkage=shrinkage, &
            return_vectors=return_vectors)
    end function varp_expected_band_multitaper_canonical_coherence

end module scifort_signal_varp_multitaper_target
