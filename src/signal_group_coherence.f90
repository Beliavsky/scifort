! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Conditional group coherence of Hermitian cross-spectral matrices.
! This is a SciFort extension; no matching scipy.signal entry point exists.
! The statistic 1-Lambda is the determinant (Wilks) summary of all canonical
! coherences between groups A and B conditional on a disjoint group C:
!   Lambda = det(S_{ABC}) det(S_C) / (det(S_{AC}) det(S_{BC})).
! A singleton A, with B the other channels, reduces to multiple coherence.
! A and B singleton, with C the remaining channels, reduces to partial coherence.
! Determinants are evaluated by Hermitian Cholesky factorizations in log space.
module scifort_signal_group_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_matrix, &
        partial_coherence_tf_matrix, partial_coherence_result, partial_coherence_tf_result, &
        partial_coherence_success, partial_coherence_invalid_input, partial_coherence_singular
    use scifort_signal_band_coherence, only : integrate_spectral_bands, &
        integrate_tf_spectral_bands
    implicit none
    private

    type, public :: group_coherence_result
        real(dp), allocatable :: frequencies(:) !! Fourier frequencies or integrated band centers
        real(dp), allocatable :: coherence(:) !! determinant group coherence, in [0,1]
        real(dp), allocatable :: log_information(:) !! -log Wilks' lambda, nonnegative
        real(dp), allocatable :: shrinkage_used(:) !! per-cell diagonal-target shrinkage
        logical, allocatable :: valid(:) !! true for positive definite cells
        integer :: status = partial_coherence_invalid_input !! success/invalid/singular
    end type group_coherence_result

    type, public :: group_coherence_tf_result
        real(dp), allocatable :: frequencies(:) !! frequency/band coordinates
        real(dp), allocatable :: times(:) !! segment center times
        real(dp), allocatable :: coherence(:,:) !! frequency/time determinant coherence
        real(dp), allocatable :: log_information(:,:) !! negative log Wilks' lambda
        real(dp), allocatable :: shrinkage_used(:,:) !! per-cell shrinkage
        logical, allocatable :: valid(:,:) !! invertible cells
        integer :: status = partial_coherence_invalid_input !! success/invalid/singular
    end type group_coherence_tf_result

    interface multitaper_group_coherence
        module procedure multitaper_group_coherence_real
        module procedure multitaper_group_coherence_complex
    end interface
    interface multitaper_tf_group_coherence
        module procedure multitaper_tf_group_coherence_real
        module procedure multitaper_tf_group_coherence_complex
    end interface
    interface multitaper_band_group_coherence
        module procedure multitaper_band_group_coherence_real
        module procedure multitaper_band_group_coherence_complex
    end interface
    interface multitaper_tf_band_group_coherence
        module procedure multitaper_tf_band_group_coherence_real
        module procedure multitaper_tf_band_group_coherence_complex
    end interface

    public :: group_coherence_matrix, group_coherence_tf_matrix
    public :: band_group_coherence_matrix, band_group_coherence_tf_matrix
    public :: multitaper_group_coherence, multitaper_tf_group_coherence
    public :: multitaper_band_group_coherence, multitaper_tf_band_group_coherence

contains

    pure logical function valid_groups(n,a,b,c) result(ok)
        integer, intent(in) :: n !! total number of spectral channels
        integer, intent(in) :: a(:) !! first channel group, one-based
        integer, intent(in) :: b(:) !! second channel group, one-based
        integer, intent(in) :: c(:) !! conditioning channels, disjoint from groups
        logical :: used(n)
        integer :: i, channel

        ok = .false.
        if (size(a) == 0 .or. size(b) == 0) return
        if (size(a)+size(b)+size(c) > n) return
        used = .false.
        do i = 1, size(a)+size(b)+size(c)
            if (i <= size(a)) then
                channel = a(i)
            else if (i <= size(a)+size(b)) then
                channel = b(i-size(a))
            else
                channel = c(i-size(a)-size(b))
            end if
            if (channel < 1 .or. channel > n) return
            if (used(channel)) return
            used(channel) = .true.
        end do
        ok = .true.
    end function valid_groups

    ! Log determinant of a principal Hermitian positive-definite submatrix.
    ! The empty conditioning block has determinant one, logdet zero.
    pure subroutine logdet_submatrix(s,indices,logdet,ok)
        complex(dp), intent(in) :: s(:,:) !! regularized Hermitian CSD
        integer, intent(in) :: indices(:) !! one-based subset of channel indices
        real(dp), intent(out) :: logdet !! natural logarithm of determinant
        logical, intent(out) :: ok !! true if Cholesky succeeded
        complex(dp) :: l(size(indices),size(indices)), term
        real(dp) :: diagonal, scale, tolerance
        integer :: i, j, k, n

        n = size(indices)
        logdet = 0.0_dp
        ok = .true.
        if (n == 0) return
        l = (0.0_dp,0.0_dp)
        scale = 0.0_dp
        do i = 1, n
            scale = max(scale,abs(s(indices(i),indices(i))))
        end do
        if (scale <= 0.0_dp) then
            ok = .false.
            return
        end if
        tolerance = 64.0_dp*epsilon(1.0_dp)*scale
        do i = 1, n
            do j = 1, i-1
                term = s(indices(i),indices(j))
                do k = 1, j-1
                    term = term - l(i,k)*conjg(l(j,k))
                end do
                l(i,j) = term/real(l(j,j),dp)
            end do
            diagonal = real(s(indices(i),indices(i)),dp)
            do k = 1, i-1
                diagonal = diagonal - abs(l(i,k))**2
            end do
            if (.not. ieee_is_finite(diagonal)) then
                ok = .false.
                return
            end if
            if (diagonal <= tolerance) then
                ok = .false.
                return
            end if
            l(i,i) = cmplx(sqrt(diagonal),0.0_dp,dp)
            logdet = logdet + log(diagonal/scale)
        end do
        logdet = logdet + real(n,dp)*log(scale)
    end subroutine logdet_submatrix

    ! Work from the same regularized spectrum used for partial coherence:
    ! off-diagonals are shrunk toward zero, then trace-scaled ridge is added.
    pure subroutine coherence_cell(s,a,b,c,ridge,alpha,coherence,information,ok)
        complex(dp), intent(in) :: s(:,:) !! unregularized channel spectral matrix
        integer, intent(in) :: a(:) !! first group
        integer, intent(in) :: b(:) !! second group
        integer, intent(in) :: c(:) !! conditioning group
        real(dp), intent(in) :: ridge !! relative trace loading
        real(dp), intent(in) :: alpha !! diagonal-target shrinkage coefficient
        real(dp), intent(out) :: coherence !! determinant summary, in [0,1]
        real(dp), intent(out) :: information !! minus log determinant ratio
        logical, intent(out) :: ok !! valid positive-definite subset
        complex(dp) :: regularized(size(s,1),size(s,2))
        integer :: i, n
        real(dp) :: logall, loga, logb, logc, logratio, mean_power, tol
        logical :: good

        ok = .false.
        coherence = 0.0_dp
        information = 0.0_dp
        n = size(s,1)
        mean_power = 0.0_dp
        do i = 1, n
            mean_power = mean_power + max(0.0_dp,real(s(i,i),dp))/real(n,dp)
        end do
        regularized = 0.5_dp*(s+transpose(conjg(s)))
        regularized = (1.0_dp-alpha)*regularized
        do i = 1, n
            regularized(i,i) = cmplx(real(s(i,i),dp)+ridge*mean_power,0.0_dp,dp)
        end do
        call logdet_submatrix(regularized,[a,b,c],logall,good)
        if (.not. good) return
        call logdet_submatrix(regularized,[a,c],loga,good)
        if (.not. good) return
        call logdet_submatrix(regularized,[b,c],logb,good)
        if (.not. good) return
        call logdet_submatrix(regularized,c,logc,good)
        if (.not. good) return
        logratio = (logall-loga)-(logb-logc)
        if (.not. ieee_is_finite(logratio)) return
        tol = 5.0e-9_dp*real(n,dp)
        if (logratio > tol) return
        ! A positive logratio of round-off magnitude corresponds to zero
        ! group coherence; do not report a negative value.
        information = max(0.0_dp,-logratio)
        ! Stable 1-exp(-information) for nearly independent groups.
        if (information < 1.0e-4_dp) then
            coherence = information*(1.0_dp-information*(0.5_dp- &
                information*(1.0_dp/6.0_dp-information/24.0_dp)))
        else
            coherence = 1.0_dp-exp(-information)
        end if
        coherence = min(1.0_dp,max(0.0_dp,coherence))
        ok = .true.
    end subroutine coherence_cell

    function group_coherence_matrix(spectrum,group_a,group_b,conditioned_on,ridge,shrinkage, &
            auto_shrinkage,effective_tapers,cv_shrinkage) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! CSD(frequency,channel,channel)
        integer, intent(in) :: group_a(:) !! one-based indices of first channel block
        integer, intent(in) :: group_b(:) !! distinct one-based indices of second channel block
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning indices
        real(dp), intent(in), optional :: ridge !! relative trace loading, nonnegative
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage, [0,1]
        logical, intent(in), optional :: auto_shrinkage !! use Gaussian-Wishart plug-in shrinkage
        real(dp), intent(in), optional :: effective_tapers !! effective independent tapered spectra
        logical, intent(in), optional :: cv_shrinkage !! use precomputed cross-validated shrinkage
        type(group_coherence_result) :: out
        type(partial_coherence_result) :: reference
        integer, allocatable :: c(:)
        integer :: nc, nf, f
        real(dp) :: lambda, n_eff, nan
        logical :: automatic, ok

        allocate(out%frequencies(0),out%coherence(0),out%log_information(0), &
            out%shrinkage_used(0),out%valid(0))
        if (spectrum%status /= spectral_status_success) return
        if (.not. allocated(spectrum%power)) return
        nc = size(spectrum%power,2)
        allocate(c(0))
        if (present(conditioned_on)) c = conditioned_on
        if (.not. valid_groups(nc,group_a,group_b,c)) return
        automatic = .false.
        if (present(auto_shrinkage)) automatic = auto_shrinkage
        if (automatic) then
            n_eff = spectrum%effective_tapers
            if (present(effective_tapers)) n_eff = effective_tapers
            reference = partial_coherence_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
                auto_shrinkage=.true.,effective_tapers=n_eff,cv_shrinkage=cv_shrinkage)
        else
            reference = partial_coherence_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
                auto_shrinkage=.false.,effective_tapers=effective_tapers,cv_shrinkage=cv_shrinkage)
        end if
        if (reference%status == partial_coherence_invalid_input) return
        nf = size(reference%frequencies)
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        deallocate(out%frequencies,out%coherence,out%log_information,out%shrinkage_used,out%valid)
        allocate(out%frequencies(nf),out%coherence(nf),out%log_information(nf), &
            out%shrinkage_used(nf),out%valid(nf))
        out%frequencies = reference%frequencies
        out%shrinkage_used = reference%shrinkage_used
        out%valid = reference%valid
        out%coherence = nan
        out%log_information = nan
        out%status = reference%status
        lambda = 0.0_dp
        if (present(ridge)) lambda = ridge
        do f = 1, nf
            if (.not. out%valid(f)) cycle
            call coherence_cell(spectrum%power(f,:,:),group_a,group_b,c,lambda, &
                out%shrinkage_used(f),out%coherence(f),out%log_information(f),ok)
            if (ok) cycle
            out%valid(f) = .false.
            out%coherence(f) = nan
            out%log_information(f) = nan
            out%status = partial_coherence_singular
        end do
    end function group_coherence_matrix

    function group_coherence_tf_matrix(spectrum,group_a,group_b,conditioned_on,ridge,shrinkage, &
            auto_shrinkage,effective_tapers,cv_shrinkage) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum !! CSD(f,t,i,j)
        integer, intent(in) :: group_a(:) !! first channel block
        integer, intent(in) :: group_b(:) !! second disjoint channel block
        integer, intent(in), optional :: conditioned_on(:) !! third disjoint channel block
        real(dp), intent(in), optional :: ridge !! relative trace loading, >= 0
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage, [0,1]
        logical, intent(in), optional :: auto_shrinkage !! Gaussian-Wishart plug-in shrinkage
        real(dp), intent(in), optional :: effective_tapers !! effective count of fixed-weight tapers
        logical, intent(in), optional :: cv_shrinkage !! use precomputed LOO shrinkage
        type(group_coherence_tf_result) :: out
        type(partial_coherence_tf_result) :: reference
        integer, allocatable :: c(:)
        integer :: nc, nf, nt, f, t
        real(dp) :: lambda, n_eff, nan
        logical :: automatic, ok

        allocate(out%frequencies(0),out%times(0),out%coherence(0,0),out%log_information(0,0), &
            out%shrinkage_used(0,0),out%valid(0,0))
        if (spectrum%status /= spectral_status_success) return
        if (.not. allocated(spectrum%power)) return
        nc = size(spectrum%power,3)
        allocate(c(0))
        if (present(conditioned_on)) c = conditioned_on
        if (.not. valid_groups(nc,group_a,group_b,c)) return
        automatic = .false.
        if (present(auto_shrinkage)) automatic = auto_shrinkage
        if (automatic) then
            n_eff = spectrum%effective_tapers
            if (present(effective_tapers)) n_eff = effective_tapers
            reference = partial_coherence_tf_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
                auto_shrinkage=.true.,effective_tapers=n_eff,cv_shrinkage=cv_shrinkage)
        else
            reference = partial_coherence_tf_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
                auto_shrinkage=.false.,effective_tapers=effective_tapers,cv_shrinkage=cv_shrinkage)
        end if
        if (reference%status == partial_coherence_invalid_input) return
        nf = size(reference%frequencies)
        nt = size(reference%times)
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        deallocate(out%frequencies,out%times,out%coherence,out%log_information,out%shrinkage_used,out%valid)
        allocate(out%frequencies(nf),out%times(nt),out%coherence(nf,nt),out%log_information(nf,nt), &
            out%shrinkage_used(nf,nt),out%valid(nf,nt))
        out%frequencies = reference%frequencies
        out%times = reference%times
        out%shrinkage_used = reference%shrinkage_used
        out%valid = reference%valid
        out%coherence = nan
        out%log_information = nan
        out%status = reference%status
        lambda = 0.0_dp
        if (present(ridge)) lambda = ridge
        do t = 1, nt
            do f = 1, nf
                if (.not. out%valid(f,t)) cycle
                call coherence_cell(spectrum%power(f,t,:,:),group_a,group_b,c,lambda, &
                    out%shrinkage_used(f,t),out%coherence(f,t),out%log_information(f,t),ok)
                if (ok) cycle
                out%valid(f,t) = .false.
                out%coherence(f,t) = nan
                out%log_information(f,t) = nan
                out%status = partial_coherence_singular
            end do
        end do
    end function group_coherence_tf_matrix

    function band_group_coherence_matrix(spectrum,bands,group_a,group_b,conditioned_on,average, &
            ridge,shrinkage) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! frequency-resolved CSD
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband), in Hz
        integer, intent(in) :: group_a(:) !! first group
        integer, intent(in) :: group_b(:) !! second group
        integer, intent(in), optional :: conditioned_on(:) !! additional conditioning group
        logical, intent(in), optional :: average !! return mean spectral density in band
        real(dp), intent(in), optional :: ridge !! relative trace loading after integration
        real(dp), intent(in), optional :: shrinkage !! diagonal shrinkage after integration
        type(group_coherence_result) :: out
        type(multitaper_matrix_result) :: integrated

        integrated = integrate_spectral_bands(spectrum,bands,average)
        out = group_coherence_matrix(integrated,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage)
    end function band_group_coherence_matrix

    function band_group_coherence_tf_matrix(spectrum,bands,group_a,group_b,conditioned_on,average, &
            ridge,shrinkage) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum !! CSD(f,t,i,j)
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband), in Hz
        integer, intent(in) :: group_a(:) !! first group
        integer, intent(in) :: group_b(:) !! second group
        integer, intent(in), optional :: conditioned_on(:) !! disjoint conditioning group
        logical, intent(in), optional :: average !! return mean spectral density in band
        real(dp), intent(in), optional :: ridge !! relative trace loading after integration
        real(dp), intent(in), optional :: shrinkage !! diagonal shrinkage after integration
        type(group_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: integrated

        integrated = integrate_tf_spectral_bands(spectrum,bands,average)
        out = group_coherence_tf_matrix(integrated,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage)
    end function band_group_coherence_tf_matrix

    function multitaper_group_coherence_real(x,nw,group_a,group_b,conditioned_on,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage,auto_shrinkage,cv_shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! real observations (sample,channel) by default
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: group_a(:) !! first group indices
        integer, intent(in) :: group_b(:) !! second group indices
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean before tapering
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! one-based sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative trace ridge
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        logical, intent(in), optional :: auto_shrinkage !! fixed-taper automatic shrinkage
        logical, intent(in), optional :: cv_shrinkage !! LOO shrinkage including adaptive weights
        type(group_coherence_result) :: out
        type(multitaper_matrix_result) :: spec

        spec = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,cv_shrinkage)
        out = group_coherence_matrix(spec,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,auto_shrinkage=auto_shrinkage,cv_shrinkage=cv_shrinkage)
    end function multitaper_group_coherence_real

    function multitaper_group_coherence_complex(x,nw,group_a,group_b,conditioned_on,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage,auto_shrinkage,cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex observations
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: group_a(:) !! first group indices
        integer, intent(in) :: group_b(:) !! second group indices
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean before tapering
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! one-based sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative trace ridge
        real(dp), intent(in), optional :: shrinkage !! fixed diagonal-target shrinkage
        logical, intent(in), optional :: auto_shrinkage !! fixed-taper automatic shrinkage
        logical, intent(in), optional :: cv_shrinkage !! LOO shrinkage including adaptive weights
        type(group_coherence_result) :: out
        type(multitaper_matrix_result) :: spec

        spec = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,cv_shrinkage)
        out = group_coherence_matrix(spec,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,auto_shrinkage=auto_shrinkage,cv_shrinkage=cv_shrinkage)
    end function multitaper_group_coherence_complex

    function multitaper_tf_group_coherence_real(x,nw,nperseg,group_a,group_b,conditioned_on, &
            noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            auto_shrinkage,cv_shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! real samples by channels
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples per window
        integer, intent(in) :: group_a(:) !! first channel group
        integer, intent(in) :: group_b(:) !! second channel group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! Fourier length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis (one-based)
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative loading
        real(dp), intent(in), optional :: shrinkage !! fixed shrinkage
        logical, intent(in), optional :: auto_shrinkage !! fixed-taper automatic shrinkage
        logical, intent(in), optional :: cv_shrinkage !! adaptive-capable predictive shrinkage
        type(group_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spec

        spec = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode,cv_shrinkage)
        out = group_coherence_tf_matrix(spec,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,auto_shrinkage=auto_shrinkage,cv_shrinkage=cv_shrinkage)
    end function multitaper_tf_group_coherence_real

    function multitaper_tf_group_coherence_complex(x,nw,nperseg,group_a,group_b,conditioned_on, &
            noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            auto_shrinkage,cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex samples by channels
        real(dp), intent(in) :: nw !! DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! samples per window
        integer, intent(in) :: group_a(:) !! first channel group
        integer, intent(in) :: group_b(:) !! second channel group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning channels
        integer, intent(in), optional :: noverlap !! overlapped samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! Fourier length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis (one-based)
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative loading
        real(dp), intent(in), optional :: shrinkage !! fixed shrinkage
        logical, intent(in), optional :: auto_shrinkage !! fixed-taper automatic shrinkage
        logical, intent(in), optional :: cv_shrinkage !! predictive shrinkage
        type(group_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spec

        spec = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode,cv_shrinkage)
        out = group_coherence_tf_matrix(spec,group_a,group_b,conditioned_on, &
            ridge=ridge,shrinkage=shrinkage,auto_shrinkage=auto_shrinkage,cv_shrinkage=cv_shrinkage)
    end function multitaper_tf_group_coherence_complex

    function multitaper_band_group_coherence_real(x,nw,bands,group_a,group_b,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! real signals
        real(dp), intent(in) :: nw !! time-bandwidth product
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband)
        integer, intent(in) :: group_a(:) !! first channel group
        integer, intent(in) :: group_b(:) !! second channel group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! Fourier length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: average !! mean PSD per band
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! shrinkage after integration
        type(group_coherence_result) :: out
        type(multitaper_matrix_result) :: spec

        spec = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_group_coherence_matrix(spec,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage)
    end function multitaper_band_group_coherence_real

    function multitaper_band_group_coherence_complex(x,nw,bands,group_a,group_b,conditioned_on, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex signals
        real(dp), intent(in) :: nw !! time-bandwidth product
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband)
        integer, intent(in) :: group_a(:) !! first channel group
        integer, intent(in) :: group_b(:) !! second channel group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! Fourier length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: average !! mean PSD per band
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! shrinkage after integration
        type(group_coherence_result) :: out
        type(multitaper_matrix_result) :: spec

        spec = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode)
        out = band_group_coherence_matrix(spec,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage)
    end function multitaper_band_group_coherence_complex

    function multitaper_tf_band_group_coherence_real(x,nw,nperseg,bands,group_a,group_b,conditioned_on, &
            noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! real signals
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! window sample count
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband)
        integer, intent(in) :: group_a(:) !! first group
        integer, intent(in) :: group_b(:) !! second group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: noverlap !! sample overlap
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: average !! mean band spectral density
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! shrinkage after integration
        type(group_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spec

        spec = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode)
        out = band_group_coherence_tf_matrix(spec,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage)
    end function multitaper_tf_band_group_coherence_real

    function multitaper_tf_band_group_coherence_complex(x,nw,nperseg,bands,group_a,group_b,conditioned_on, &
            noverlap,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,average,ridge,shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex signals
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! window sample count
        real(dp), intent(in) :: bands(:,:) !! band edges (2,nband)
        integer, intent(in) :: group_a(:) !! first group
        integer, intent(in) :: group_b(:) !! second group
        integer, intent(in), optional :: conditioned_on(:) !! conditioning group
        integer, intent(in), optional :: noverlap !! sample overlap
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! sample axis
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: average !! mean band spectral density
        real(dp), intent(in), optional :: ridge !! relative loading after integration
        real(dp), intent(in), optional :: shrinkage !! shrinkage after integration
        type(group_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spec

        spec = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode)
        out = band_group_coherence_tf_matrix(spec,bands,group_a,group_b,conditioned_on, &
            average,ridge,shrinkage)
    end function multitaper_tf_band_group_coherence_complex

end module scifort_signal_group_coherence
