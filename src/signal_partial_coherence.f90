! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Conditional (partial) coherence and complex coherency from spectral matrices.
! This SciFort extension uses a Hermitian Cholesky factorization at each FFT
! frequency, avoiding explicitly expanded determinant/cofactor formulas.
! Optional diagonal-target shrinkage stabilizes poorly conditioned spectra.
module scifort_signal_partial_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success
    use scifort_signal_multitaper, only : multitaper_matrix_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix_spectrogram
    implicit none
    private

    integer, parameter, public :: partial_coherence_success = 0
    integer, parameter, public :: partial_coherence_invalid_input = 1
    integer, parameter, public :: partial_coherence_singular = 2

    type, public :: partial_coherence_result
        real(dp), allocatable :: frequencies(:) !! frequency grid, copied from input
        real(dp), allocatable :: coherence(:,:,:) !! squared partial coherence (frequency,channel,channel)
        complex(dp), allocatable :: coherency(:,:,:) !! optional complex partial coherency, when requested
        complex(dp), allocatable :: precision(:,:,:) !! inverse regularized spectral matrices
        logical, allocatable :: valid(:) !! invertible frequency bins
        real(dp), allocatable :: shrinkage_used(:) !! diagonal shrinkage chosen at each frequency
        integer :: status = partial_coherence_invalid_input !! success, invalid input, or singular bins
    end type partial_coherence_result

    ! Time-resolved analogue. Every time-frequency matrix is inverted
    ! independently; the number of tapers limits its unregularized rank.
    type, public :: partial_coherence_tf_result
        real(dp), allocatable :: frequencies(:) !! Fourier frequencies
        real(dp), allocatable :: times(:) !! segment-center times
        real(dp), allocatable :: coherence(:,:,:,:) !! frequency, time, channel, channel
        complex(dp), allocatable :: coherency(:,:,:,:) !! optional complex partial coherency (f,t,i,j)
        complex(dp), allocatable :: precision(:,:,:,:) !! regularized inverse CSD matrices
        logical, allocatable :: valid(:,:) !! frequency, time invertibility mask
        real(dp), allocatable :: shrinkage_used(:,:) !! chosen shrinkage at each frequency and time
        integer :: status = partial_coherence_invalid_input !! success, invalid, or singular
    end type partial_coherence_tf_result

    interface multitaper_tf_partial_coherence
        module procedure multitaper_tf_partial_coherence_real
        module procedure multitaper_tf_partial_coherence_complex
    end interface multitaper_tf_partial_coherence

    interface multitaper_partial_coherence
        module procedure multitaper_partial_coherence_real
        module procedure multitaper_partial_coherence_complex
    end interface multitaper_partial_coherence

    public :: partial_coherence_matrix, multitaper_partial_coherence
    public :: partial_coherence_tf_matrix, multitaper_tf_partial_coherence

contains

    function partial_coherence_matrix(spectrum, ridge, shrinkage, return_coherency, &
            auto_shrinkage, effective_tapers, cv_shrinkage) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! Hermitian cross-spectrum matrices
        real(dp), intent(in), optional :: ridge !! nonnegative relative diagonal loading, default zero
        real(dp), intent(in), optional :: shrinkage !! diagonal-target shrinkage in [0,1], default zero
        logical, intent(in), optional :: return_coherency !! compute complex partial coherency if true
        logical, intent(in), optional :: auto_shrinkage !! estimate per-frequency diagonal shrinkage
        real(dp), intent(in), optional :: effective_tapers !! effective independent complex samples >= 1
        logical, intent(in), optional :: cv_shrinkage !! apply LOO coefficients stored in spectrum
        type(partial_coherence_result) :: out

        complex(dp), allocatable :: a(:,:), l(:,:), inverse(:,:), work(:)
        real(dp) :: lambda, alpha, scale, norm_a, tol, pivot, value, nan, n_effective, local_alpha
        complex(dp) :: term
        integer :: nf, nc, f, i, j, k
        logical :: failed, with_coherency, automatic, cross_validate

        allocate(out%frequencies(0), out%coherence(0,0,0), out%coherency(0,0,0), &
            out%precision(0,0,0), out%valid(0), out%shrinkage_used(0))
        if (spectrum%status /= spectral_status_success) return
        if (.not. allocated(spectrum%frequencies)) return
        if (.not. allocated(spectrum%power)) return
        nf = size(spectrum%power,1)
        nc = size(spectrum%power,2)
        if (nf < 1 .or. nc < 2 .or. size(spectrum%power,3) /= nc) return
        if (size(spectrum%frequencies) /= nf) return
        if (.not. all(ieee_is_finite(spectrum%frequencies))) return
        if (.not. all(ieee_is_finite(real(spectrum%power,dp)))) return
        if (.not. all(ieee_is_finite(aimag(spectrum%power)))) return
        lambda = 0.0_dp
        if (present(ridge)) lambda = ridge
        if (.not. ieee_is_finite(lambda)) return
        if (lambda < 0.0_dp) return
        alpha = 0.0_dp
        if (present(shrinkage)) alpha = shrinkage
        if (.not. ieee_is_finite(alpha)) return
        if (alpha < 0.0_dp .or. alpha > 1.0_dp) return
        automatic = .false.
        if (present(auto_shrinkage)) automatic = auto_shrinkage
        n_effective = 0.0_dp
        if (automatic) then
            if (present(shrinkage)) return
            if (.not. present(effective_tapers)) return
            n_effective = effective_tapers
            if (.not. ieee_is_finite(n_effective)) return
            if (n_effective < 1.0_dp) return
        else
            if (present(effective_tapers)) return
        end if
        cross_validate = .false.
        if (present(cv_shrinkage)) cross_validate = cv_shrinkage
        if (cross_validate) then
            if (automatic .or. present(shrinkage) .or. present(effective_tapers)) return
            if (.not. allocated(spectrum%cv_shrinkage)) return
            if (size(spectrum%cv_shrinkage) /= nf) return
            if (.not. all(ieee_is_finite(spectrum%cv_shrinkage))) return
            if (any(spectrum%cv_shrinkage < 0.0_dp .or. spectrum%cv_shrinkage > 1.0_dp)) return
        end if
        with_coherency = .false.
        if (present(return_coherency)) with_coherency = return_coherency

        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        deallocate(out%frequencies, out%coherence, out%precision, out%valid, out%shrinkage_used)
        allocate(out%frequencies(nf), out%coherence(nf,nc,nc), &
            out%precision(nf,nc,nc), out%valid(nf), out%shrinkage_used(nf))
        if (with_coherency) then
            deallocate(out%coherency)
            allocate(out%coherency(nf,nc,nc))
            out%coherency = cmplx(nan,nan,dp)
        end if
        out%frequencies = spectrum%frequencies
        out%coherence = nan
        out%precision = cmplx(nan,nan,dp)
        out%valid = .false.
        out%shrinkage_used = nan
        out%status = partial_coherence_success
        allocate(a(nc,nc),l(nc,nc),inverse(nc,nc),work(nc))

        do f = 1, nf
            a = spectrum%power(f,:,:)
            norm_a = maxval(abs(a))
            if (norm_a <= 0.0_dp) then
                out%status = partial_coherence_singular
                cycle
            end if
            tol = 64.0_dp*epsilon(1.0_dp)*norm_a
            ! The supplied matrix must actually represent a Hermitian PSD.
            ! Refuse non-Hermitian inputs rather than silently misinterpreting them.
            if (maxval(abs(a-transpose(conjg(a)))) > 16.0_dp*tol) then
                out%status = partial_coherence_invalid_input
                return
            end if
            a = 0.5_dp*(a+transpose(conjg(a)))
            scale = 0.0_dp
            do i = 1, nc
                if (real(a(i,i),dp) < -tol) then
                    out%status = partial_coherence_invalid_input
                    return
                end if
                scale = scale + max(0.0_dp,real(a(i,i),dp))
            end do
            scale = scale/real(nc,dp)
            if (scale <= 0.0_dp) then
                out%status = partial_coherence_singular
                cycle
            end if
            ! Gaussian-Wishart plug-in shrinkage is evaluated on the ORIGINAL
            ! Hermitian matrix (before adding ridge); it uses scale-free ratios.
            local_alpha = alpha
            if (automatic) local_alpha = automatic_diagonal_shrinkage(a, n_effective)
            if (cross_validate) local_alpha = spectrum%cv_shrinkage(f)
            out%shrinkage_used(f) = local_alpha
            ! Diagonal-target shrinkage retains each channel power, reducing
            ! cross spectra by (1-alpha). Unlike trace ridge, it commutes
            ! with independent rescaling of individual channels.
            do j = 1, nc
                do i = 1, nc
                    if (i /= j) a(i,j) = (1.0_dp-local_alpha)*a(i,j)
                end do
            end do
            ! Trace-based ridge loading is dimensionless in common power units.
            do i = 1, nc
                a(i,i) = a(i,i) + cmplx(lambda*scale,0.0_dp,dp)
            end do
            tol = 64.0_dp*epsilon(1.0_dp)*maxval(abs(a))
            l = (0.0_dp,0.0_dp)
            failed = .false.
            do i = 1, nc
                do j = 1, i-1
                    term = a(i,j)
                    do k = 1, j-1
                        term = term - l(i,k)*conjg(l(j,k))
                    end do
                    l(i,j) = term/real(l(j,j),dp)
                end do
                pivot = real(a(i,i),dp)
                do k = 1, i-1
                    pivot = pivot - abs(l(i,k))**2
                end do
                if (.not. ieee_is_finite(pivot)) then
                    failed = .true.
                    exit
                end if
                if (pivot <= tol) then
                    failed = .true.
                    exit
                end if
                l(i,i) = cmplx(sqrt(pivot),0.0_dp,dp)
            end do
            if (failed) then
                out%status = partial_coherence_singular
                cycle
            end if

            ! Solve L*L^H * inverse = I by forward/back substitution.
            inverse = (0.0_dp,0.0_dp)
            do j = 1, nc
                do i = 1, nc
                    work(i) = (0.0_dp,0.0_dp)
                    if (i == j) work(i) = (1.0_dp,0.0_dp)
                    do k = 1, i-1
                        work(i) = work(i) - l(i,k)*work(k)
                    end do
                    work(i) = work(i)/real(l(i,i),dp)
                end do
                do i = nc, 1, -1
                    inverse(i,j) = work(i)
                    do k = i+1, nc
                        inverse(i,j) = inverse(i,j) - conjg(l(k,i))*inverse(k,j)
                    end do
                    inverse(i,j) = inverse(i,j)/real(l(i,i),dp)
                end do
            end do
            inverse = 0.5_dp*(inverse+transpose(conjg(inverse)))
            failed = .false.
            do i = 1, nc
                if (.not. ieee_is_finite(real(inverse(i,i),dp))) failed = .true.
                if (real(inverse(i,i),dp) <= 0.0_dp) failed = .true.
            end do
            if (failed) then
                out%status = partial_coherence_singular
                cycle
            end if
            out%precision(f,:,:) = inverse
            do j = 1, nc
                do i = 1, nc
                    ! Normalize before squaring to avoid overflow for
                    ! small spectral-density units or large precision values.
                    value = abs(inverse(i,j))/sqrt(real(inverse(i,i),dp))
                    value = (value/sqrt(real(inverse(j,j),dp)))**2
                    if (.not. ieee_is_finite(value)) then
                        failed = .true.
                        exit
                    end if
                    out%coherence(f,i,j) = min(1.0_dp,max(0.0_dp,value))
                    if (with_coherency) then
                        if (i == j) then
                            out%coherency(f,i,j) = (1.0_dp,0.0_dp)
                        else
                            term = -inverse(i,j)/sqrt(real(inverse(i,i),dp))
                            term = term/sqrt(real(inverse(j,j),dp))
                            if (.not. ieee_is_finite(real(term,dp)) .or. &
                                .not. ieee_is_finite(aimag(term))) then
                                failed = .true.
                                exit
                            end if
                            out%coherency(f,i,j) = term
                        end if
                    end if
                end do
                if (failed) exit
            end do
            if (failed) then
                out%coherence(f,:,:) = nan
                out%precision(f,:,:) = cmplx(nan,nan,dp)
                if (with_coherency) out%coherency(f,:,:) = cmplx(nan,nan,dp)
                out%status = partial_coherence_singular
                cycle
            end if
            do i = 1, nc
                out%coherence(f,i,i) = 1.0_dp
            end do
            out%valid(f) = .true.
        end do
    end function partial_coherence_matrix

    ! Gaussian-Wishart diagonal-target estimate. For n independent proper
    ! complex Gaussian tapered Fourier vectors, S is a complex Wishart sample
    ! covariance. Its off-diagonal sampling variance is Sigma_ii*Sigma_jj/n.
    ! The moment identities E[Sii*Sjj] = Sigma_ii*Sigma_jj + |Sigma_ij|^2/n
    ! and E[|Sij|^2] = |Sigma_ij|^2 + Sigma_ii*Sigma_jj/n yield
    !     alpha = clip((n*D-O) / ((n*n-1)*O), 0, 1),
    ! with D=sum(i/=j) Sii*Sjj and O=sum(i/=j) |Sij|^2.
    ! This is a plug-in, not exact oracle or classical Ledoit-Wolf shrinkage.
    ! For tapered or colored data, effective_tapers is an approximation;
    ! adjacent tapers need not be independent or proper-complex at DC/Nyquist.
    pure function automatic_diagonal_shrinkage(s, n) result(alpha)
        complex(dp), intent(in) :: s(:,:) !! Hermitian PSD matrix
        real(dp), intent(in) :: n !! effective number of independent complex observations >= 1
        real(dp) :: alpha, diagonal_product, offdiagonal_power, scale, d, o
        integer :: i, j, nc

        nc = size(s,1)
        if (n <= 1.0_dp+epsilon(n)) then
            alpha = 1.0_dp
            return
        end if
        scale = maxval(abs(s))
        diagonal_product = 0.0_dp
        offdiagonal_power = 0.0_dp
        do j = 1, nc
            do i = 1, nc
                if (i == j) cycle
                d = max(0.0_dp,real(s(i,i),dp)/scale)*max(0.0_dp,real(s(j,j),dp)/scale)
                o = abs(s(i,j)/scale)**2
                diagonal_product = diagonal_product+d
                offdiagonal_power = offdiagonal_power+o
            end do
        end do
        if (offdiagonal_power == 0.0_dp) then
            alpha = 0.0_dp
        else
            alpha = max(0.0_dp,min(1.0_dp, &
                (n*diagonal_product-offdiagonal_power)/ &
                ((n-1.0_dp)*(n+1.0_dp)*offdiagonal_power)))
        end if
    end function automatic_diagonal_shrinkage

    function multitaper_partial_coherence_real(x, nw, kmax, fs, nfft, detrend, &
            weighting, axis, detrend_mode, ridge, shrinkage, return_coherency, auto_shrinkage, cv_shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! observations, samples x channels by default
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample frequency
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! enable mean removal
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative diagonal loading
        real(dp), intent(in), optional :: shrinkage !! diagonal-target shrinkage in [0,1]
        logical, intent(in), optional :: return_coherency !! return complex coherency when true
        logical, intent(in), optional :: auto_shrinkage !! estimate diagonal shrinkage from effective taper count
        logical, intent(in), optional :: cv_shrinkage !! leave-one-taper-out predictive shrinkage
        type(partial_coherence_result) :: out
        type(multitaper_matrix_result) :: spectrum
        spectrum = multitaper_spectral_matrix(x, nw, kmax, fs, nfft, detrend, &
            weighting, axis, detrend_mode, cv_shrinkage)
        if (present(cv_shrinkage)) then
            if (cv_shrinkage) then
                out = partial_coherence_matrix(spectrum, ridge=ridge, shrinkage=shrinkage, &
                    return_coherency=return_coherency, auto_shrinkage=auto_shrinkage, cv_shrinkage=.true.)
                return
            end if
        end if
        if (present(auto_shrinkage)) then
            if (auto_shrinkage) then
                ! Eigen DPSS weights are deterministic, with equivalent
                ! weighted Wishart sample count; adaptive weights are not.
                if (present(weighting)) then
                    if (trim(weighting) == 'adaptive') then
                        out = partial_coherence_matrix(spectrum, auto_shrinkage=.true.)
                        return
                    end if
                end if
                if (present(shrinkage)) then
                    out = partial_coherence_matrix(spectrum, auto_shrinkage=.true.)
                    return
                end if
                out = partial_coherence_matrix(spectrum, ridge=ridge, return_coherency=return_coherency, &
                    auto_shrinkage=.true., effective_tapers=spectrum%effective_tapers)
                return
            end if
        end if
        out = partial_coherence_matrix(spectrum, ridge, shrinkage, return_coherency)
    end function multitaper_partial_coherence_real

    function multitaper_partial_coherence_complex(x, nw, kmax, fs, nfft, detrend, &
            weighting, axis, detrend_mode, ridge, shrinkage, return_coherency, auto_shrinkage, cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex observations, samples x channels by default
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample frequency
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! enable mean removal
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! relative diagonal loading
        real(dp), intent(in), optional :: shrinkage !! diagonal-target shrinkage in [0,1]
        logical, intent(in), optional :: return_coherency !! return complex coherency when true
        logical, intent(in), optional :: auto_shrinkage !! estimate diagonal shrinkage from effective taper count
        logical, intent(in), optional :: cv_shrinkage !! leave-one-taper-out predictive shrinkage
        type(partial_coherence_result) :: out
        type(multitaper_matrix_result) :: spectrum
        spectrum = multitaper_spectral_matrix(x, nw, kmax, fs, nfft, detrend, &
            weighting, axis, detrend_mode, cv_shrinkage)
        if (present(cv_shrinkage)) then
            if (cv_shrinkage) then
                out = partial_coherence_matrix(spectrum, ridge=ridge, shrinkage=shrinkage, &
                    return_coherency=return_coherency, auto_shrinkage=auto_shrinkage, cv_shrinkage=.true.)
                return
            end if
        end if
        if (present(auto_shrinkage)) then
            if (auto_shrinkage) then
                ! Eigen DPSS weights are deterministic, with equivalent
                ! weighted Wishart sample count; adaptive weights are not.
                if (present(weighting)) then
                    if (trim(weighting) == 'adaptive') then
                        out = partial_coherence_matrix(spectrum, auto_shrinkage=.true.)
                        return
                    end if
                end if
                if (present(shrinkage)) then
                    out = partial_coherence_matrix(spectrum, auto_shrinkage=.true.)
                    return
                end if
                out = partial_coherence_matrix(spectrum, ridge=ridge, return_coherency=return_coherency, &
                    auto_shrinkage=.true., effective_tapers=spectrum%effective_tapers)
                return
            end if
        end if
        out = partial_coherence_matrix(spectrum, ridge, shrinkage, return_coherency)
    end function multitaper_partial_coherence_complex

    ! Invert every channel spectral matrix at every time-frequency bin. Reuse
    ! the existing precision/regularization kernel rather than maintaining a
    ! second Cholesky implementation whose singularity rules might diverge.
    function partial_coherence_tf_matrix(spectrum, ridge, shrinkage, return_coherency, &
            auto_shrinkage, effective_tapers, cv_shrinkage) result(out)
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum !! CSD(f,t,i,j)
        real(dp), intent(in), optional :: ridge !! relative trace-based diagonal loading >= 0
        real(dp), intent(in), optional :: shrinkage !! diagonal-target shrinkage in [0,1]
        logical, intent(in), optional :: return_coherency !! return complex coherency when true
        logical, intent(in), optional :: auto_shrinkage !! estimate diagonal shrinkage from effective taper count
        real(dp), intent(in), optional :: effective_tapers !! independent complex observations >= 1
        logical, intent(in), optional :: cv_shrinkage !! apply LOO coefficients stored in spectrum
        type(partial_coherence_tf_result) :: out
        type(multitaper_matrix_result) :: frame
        type(partial_coherence_result) :: result_frame
        integer :: nf, nt, nc, t
        real(dp) :: nan, lambda, alpha, n_effective
        logical :: with_coherency, automatic, cross_validate

        allocate(out%frequencies(0), out%times(0), &
            out%coherence(0,0,0,0), out%coherency(0,0,0,0), &
            out%precision(0,0,0,0), out%valid(0,0), out%shrinkage_used(0,0))
        if (spectrum%status /= spectral_status_success) return
        if (.not. allocated(spectrum%frequencies)) return
        if (.not. allocated(spectrum%times)) return
        if (.not. allocated(spectrum%power)) return
        nf = size(spectrum%power, 1)
        nt = size(spectrum%power, 2)
        nc = size(spectrum%power, 3)
        if (nf < 1 .or. nt < 1 .or. nc < 2) return
        if (size(spectrum%power, 4) /= nc) return
        if (size(spectrum%frequencies) /= nf) return
        if (size(spectrum%times) /= nt) return
        if (.not. all(ieee_is_finite(spectrum%frequencies))) return
        if (.not. all(ieee_is_finite(spectrum%times))) return
        if (.not. all(ieee_is_finite(real(spectrum%power, dp)))) return
        if (.not. all(ieee_is_finite(aimag(spectrum%power)))) return
        lambda = 0.0_dp
        if (present(ridge)) lambda = ridge
        if (.not. ieee_is_finite(lambda)) return
        if (lambda < 0.0_dp) return
        alpha = 0.0_dp
        if (present(shrinkage)) alpha = shrinkage
        if (.not. ieee_is_finite(alpha)) return
        if (alpha < 0.0_dp .or. alpha > 1.0_dp) return
        automatic = .false.
        if (present(auto_shrinkage)) automatic = auto_shrinkage
        n_effective = 0.0_dp
        if (automatic) then
            if (present(shrinkage)) return
            if (.not. present(effective_tapers)) return
            n_effective = effective_tapers
            if (.not. ieee_is_finite(n_effective)) return
            if (n_effective < 1.0_dp) return
        else
            if (present(effective_tapers)) return
        end if
        cross_validate = .false.
        if (present(cv_shrinkage)) cross_validate = cv_shrinkage
        if (cross_validate) then
            if (automatic .or. present(shrinkage) .or. present(effective_tapers)) return
            if (.not. allocated(spectrum%cv_shrinkage)) return
            if (size(spectrum%cv_shrinkage,1) /= nf .or. &
                size(spectrum%cv_shrinkage,2) /= nt) return
            if (.not. all(ieee_is_finite(spectrum%cv_shrinkage))) return
            if (any(spectrum%cv_shrinkage < 0.0_dp .or. spectrum%cv_shrinkage > 1.0_dp)) return
        end if
        with_coherency = .false.
        if (present(return_coherency)) with_coherency = return_coherency

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        deallocate(out%frequencies, out%times, out%coherence, &
            out%precision, out%valid, out%shrinkage_used)
        allocate(out%frequencies(nf), out%times(nt), &
            out%coherence(nf,nt,nc,nc), out%precision(nf,nt,nc,nc), &
            out%valid(nf,nt), out%shrinkage_used(nf,nt))
        if (with_coherency) then
            deallocate(out%coherency)
            allocate(out%coherency(nf,nt,nc,nc))
            out%coherency = cmplx(nan,nan,dp)
        end if
        out%frequencies = spectrum%frequencies
        out%times = spectrum%times
        out%coherence = nan
        out%precision = cmplx(nan, nan, dp)
        out%valid = .false.
        out%shrinkage_used = nan
        out%status = partial_coherence_success

        allocate(frame%frequencies(nf), frame%power(nf,nc,nc))
        frame%status = spectral_status_success
        frame%frequencies = spectrum%frequencies
        do t = 1, nt
            frame%power = spectrum%power(:,t,:,:)
            if (cross_validate) frame%cv_shrinkage = spectrum%cv_shrinkage(:,t)
            if (cross_validate) then
                result_frame = partial_coherence_matrix(frame, ridge=lambda, &
                    cv_shrinkage=.true., return_coherency=with_coherency)
            else if (automatic) then
                result_frame = partial_coherence_matrix(frame, ridge=lambda, &
                    auto_shrinkage=.true., effective_tapers=n_effective, &
                    return_coherency=with_coherency)
            else
                result_frame = partial_coherence_matrix(frame, ridge=lambda, &
                    shrinkage=alpha, return_coherency=with_coherency)
            end if
            if (result_frame%status == partial_coherence_invalid_input) then
                out%status = partial_coherence_invalid_input
                return
            end if
            if (result_frame%status == partial_coherence_singular) &
                out%status = partial_coherence_singular
            out%coherence(:,t,:,:) = result_frame%coherence
            if (with_coherency) out%coherency(:,t,:,:) = result_frame%coherency
            out%precision(:,t,:,:) = result_frame%precision
            out%valid(:,t) = result_frame%valid
            out%shrinkage_used(:,t) = result_frame%shrinkage_used
        end do
    end function partial_coherence_tf_matrix

    function multitaper_tf_partial_coherence_real(x, nw, nperseg, noverlap, kmax, fs, &
            nfft, detrend, weighting, axis, detrend_mode, ridge, shrinkage, return_coherency, &
            auto_shrinkage, cv_shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! observed real samples, samples x channels by default
        real(dp), intent(in) :: nw !! DPSS time-half-bandwidth product
        integer, intent(in) :: nperseg !! number of samples per moving segment
        integer, intent(in), optional :: noverlap !! number of overlapping samples
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! positive sample rate
        integer, intent(in), optional :: nfft !! Fourier length >= nperseg
        logical, intent(in), optional :: detrend !! whether to remove segment means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive taper weights
        integer, intent(in), optional :: axis !! sample axis, either 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! nonnegative relative diagonal loading
        real(dp), intent(in), optional :: shrinkage !! diagonal-target shrinkage in [0,1]
        logical, intent(in), optional :: return_coherency !! return complex coherency when true
        logical, intent(in), optional :: auto_shrinkage !! estimate diagonal shrinkage from effective taper count
        logical, intent(in), optional :: cv_shrinkage !! leave-one-taper-out predictive shrinkage
        type(partial_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectrum

        spectrum = multitaper_spectral_matrix_spectrogram(x, nw, nperseg, noverlap, &
            kmax, fs, nfft, detrend, weighting, axis, detrend_mode, cv_shrinkage)
        if (present(cv_shrinkage)) then
            if (cv_shrinkage) then
                out = partial_coherence_tf_matrix(spectrum, ridge=ridge, shrinkage=shrinkage, &
                    return_coherency=return_coherency, auto_shrinkage=auto_shrinkage, cv_shrinkage=.true.)
                return
            end if
        end if
        if (present(auto_shrinkage)) then
            if (auto_shrinkage) then
                if (present(weighting)) then
                    if (trim(weighting) == 'adaptive') then
                        out = partial_coherence_tf_matrix(spectrum, auto_shrinkage=.true.)
                        return
                    end if
                end if
                if (present(shrinkage)) then
                    out = partial_coherence_tf_matrix(spectrum, auto_shrinkage=.true.)
                    return
                end if
                out = partial_coherence_tf_matrix(spectrum, ridge=ridge, return_coherency=return_coherency, &
                    auto_shrinkage=.true., effective_tapers=spectrum%effective_tapers)
                return
            end if
        end if
        out = partial_coherence_tf_matrix(spectrum, ridge, shrinkage, return_coherency)
    end function multitaper_tf_partial_coherence_real

    function multitaper_tf_partial_coherence_complex(x, nw, nperseg, noverlap, kmax, fs, &
            nfft, detrend, weighting, axis, detrend_mode, ridge, shrinkage, return_coherency, &
            auto_shrinkage, cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! observed complex samples, samples x channels by default
        real(dp), intent(in) :: nw !! DPSS time-half-bandwidth product
        integer, intent(in) :: nperseg !! number of samples per moving segment
        integer, intent(in), optional :: noverlap !! number of overlapping samples
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! positive sample rate
        integer, intent(in), optional :: nfft !! Fourier length >= nperseg
        logical, intent(in), optional :: detrend !! whether to remove segment means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive taper weights
        integer, intent(in), optional :: axis !! sample axis, either 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        real(dp), intent(in), optional :: ridge !! nonnegative relative diagonal loading
        real(dp), intent(in), optional :: shrinkage !! diagonal-target shrinkage in [0,1]
        logical, intent(in), optional :: return_coherency !! return complex coherency when true
        logical, intent(in), optional :: auto_shrinkage !! estimate diagonal shrinkage from effective taper count
        logical, intent(in), optional :: cv_shrinkage !! leave-one-taper-out predictive shrinkage
        type(partial_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectrum

        spectrum = multitaper_spectral_matrix_spectrogram(x, nw, nperseg, noverlap, &
            kmax, fs, nfft, detrend, weighting, axis, detrend_mode, cv_shrinkage)
        if (present(cv_shrinkage)) then
            if (cv_shrinkage) then
                out = partial_coherence_tf_matrix(spectrum, ridge=ridge, shrinkage=shrinkage, &
                    return_coherency=return_coherency, auto_shrinkage=auto_shrinkage, cv_shrinkage=.true.)
                return
            end if
        end if
        if (present(auto_shrinkage)) then
            if (auto_shrinkage) then
                if (present(weighting)) then
                    if (trim(weighting) == 'adaptive') then
                        out = partial_coherence_tf_matrix(spectrum, auto_shrinkage=.true.)
                        return
                    end if
                end if
                if (present(shrinkage)) then
                    out = partial_coherence_tf_matrix(spectrum, auto_shrinkage=.true.)
                    return
                end if
                out = partial_coherence_tf_matrix(spectrum, ridge=ridge, return_coherency=return_coherency, &
                    auto_shrinkage=.true., effective_tapers=spectrum%effective_tapers)
                return
            end if
        end if
        out = partial_coherence_tf_matrix(spectrum, ridge, shrinkage, return_coherency)
    end function multitaper_tf_partial_coherence_complex

end module scifort_signal_partial_coherence
