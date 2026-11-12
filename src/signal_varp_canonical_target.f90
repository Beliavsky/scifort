! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Conditional canonical coherence population targets derived from the exact
! finite-taper VAR(p) spectrum or its analytic band integrals.
! Original SciFort extension; no equivalent scipy.signal function.
module scifort_signal_varp_canonical_target
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_matrix_result
    use scifort_signal_varp_taper_spectrum, only : varp_taper_spectrum_result, &
        varp_expected_taper_spectrum, varp_expected_band_taper_spectrum, varp_spectrum_success
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, canonical_coherence_matrix
    implicit none
    private

    public :: varp_expected_canonical_coherence, varp_expected_band_canonical_coherence

contains

    function varp_expected_canonical_coherence(a,innovation_covariance,taper,frequencies, &
            group_a,group_b,conditioned_on,ridge,shrinkage,normalized,return_vectors) result(out)
        real(dp), intent(in) :: a(:,:,:) !! stable real VAR(p) coefficients (n,n,p)
        real(dp), intent(in) :: innovation_covariance(:,:) !! positive-definite innovation covariance
        real(dp), intent(in) :: taper(:) !! finite time-window weights
        real(dp), intent(in) :: frequencies(:) !! target frequencies in cycles/sample
        integer, intent(in) :: group_a(:) !! nonempty one-based channel group A
        integer, intent(in) :: group_b(:) !! nonempty disjoint one-based channel group B
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels; empty by default
        real(dp), intent(in), optional :: ridge !! nonnegative relative spectral-matrix ridge loading
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal shrinkage fraction in [0,1]
        logical, intent(in), optional :: normalized !! normalize taper energy; default false
        logical, intent(in), optional :: return_vectors !! include canonical coefficients; default false
        type(canonical_coherence_result) :: out
        type(varp_taper_spectrum_result) :: exact
        type(multitaper_matrix_result) :: s

        exact = varp_expected_taper_spectrum(a,innovation_covariance,taper,frequencies,normalized)
        if (exact%status /= varp_spectrum_success) return
        s%frequencies = exact%frequencies
        s%power = exact%power
        s%status = spectral_status_success
        out = canonical_coherence_matrix(s,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,return_vectors=return_vectors)
    end function varp_expected_canonical_coherence

    function varp_expected_band_canonical_coherence(a,innovation_covariance,taper,bands, &
            group_a,group_b,conditioned_on,average,ridge,shrinkage,normalized,return_vectors) result(out)
        real(dp), intent(in) :: a(:,:,:) !! stable real VAR(p) coefficients (n,n,p)
        real(dp), intent(in) :: innovation_covariance(:,:) !! positive-definite innovation covariance
        real(dp), intent(in) :: taper(:) !! finite time-window weights
        real(dp), intent(in) :: bands(:,:) !! exact band boundaries (2,nband), in cycles/sample
        integer, intent(in) :: group_a(:) !! nonempty one-based channel group A
        integer, intent(in) :: group_b(:) !! nonempty disjoint one-based channel group B
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning channels
        logical, intent(in), optional :: average !! average rather than integrate within bands; default false
        real(dp), intent(in), optional :: ridge !! relative spectral ridge loading, nonnegative
        real(dp), intent(in), optional :: shrinkage !! diagonal shrinkage in [0,1]
        logical, intent(in), optional :: normalized !! divide by taper energy; default false
        logical, intent(in), optional :: return_vectors !! return canonical coefficient vectors; default false
        type(canonical_coherence_result) :: out
        type(varp_taper_spectrum_result) :: exact
        type(multitaper_matrix_result) :: s

        exact = varp_expected_band_taper_spectrum(a,innovation_covariance,taper,bands, &
            average=average,normalized=normalized)
        if (exact%status /= varp_spectrum_success) return
        s%frequencies = exact%frequencies
        s%power = exact%power
        s%status = spectral_status_success
        out = canonical_coherence_matrix(s,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,return_vectors=return_vectors)
    end function varp_expected_band_canonical_coherence

end module scifort_signal_varp_canonical_target
