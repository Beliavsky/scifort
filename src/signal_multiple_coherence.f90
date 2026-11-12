! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Spectral multiple coherence: fraction of a channel's regularized power
! predictable from all other channels jointly (SciFort extension).
! Reuse the precision solver and shrinkage rules of partial coherence.
module scifort_signal_multiple_coherence
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_spectral_matrix, &
        multitaper_spectral_matrix_spectrogram
    use scifort_signal_partial_coherence, only : partial_coherence_matrix, &
        partial_coherence_tf_matrix, partial_coherence_result, partial_coherence_tf_result, &
        partial_coherence_success, partial_coherence_invalid_input, partial_coherence_singular
    implicit none
    private

    type, public :: multiple_coherence_result
        real(dp), allocatable :: frequencies(:) !! Fourier frequencies, same as input
        real(dp), allocatable :: coherence(:,:) !! squared multiple coherence (frequency, channel)
        real(dp), allocatable :: conditional_power(:,:) !! residual PSD given other channels
        real(dp), allocatable :: explained_power(:,:) !! fitted power on the regularized scale
        real(dp), allocatable :: shrinkage_used(:) !! selected diagonal shrinkage
        logical, allocatable :: valid(:) !! invertibility of each conditioned spectral matrix
        integer :: status = partial_coherence_invalid_input !! inversion status
    end type multiple_coherence_result

    type, public :: multiple_coherence_tf_result
        real(dp), allocatable :: frequencies(:) !! Fourier frequencies
        real(dp), allocatable :: times(:) !! segment-center times
        real(dp), allocatable :: coherence(:,:,:) !! multiple coherence, (frequency,time,channel)
        real(dp), allocatable :: conditional_power(:,:,:) !! residual PSD (frequency,time,channel)
        real(dp), allocatable :: explained_power(:,:,:) !! explained PSD of regularized spectrum
        real(dp), allocatable :: shrinkage_used(:,:) !! shrinkage chosen per frequency and time
        logical, allocatable :: valid(:,:) !! invertibility of each spectral matrix
        integer :: status = partial_coherence_invalid_input !! inversion status
    end type multiple_coherence_tf_result

    interface multitaper_multiple_coherence
        module procedure multitaper_multiple_coherence_real
        module procedure multitaper_multiple_coherence_complex
    end interface
    interface multitaper_tf_multiple_coherence
        module procedure multitaper_tf_multiple_coherence_real
        module procedure multitaper_tf_multiple_coherence_complex
    end interface
    public :: multiple_coherence_matrix, multiple_coherence_tf_matrix
    public :: multitaper_multiple_coherence, multitaper_tf_multiple_coherence

contains

    ! The precision-matrix identity avoids repeated regressions and remains
    ! meaningful even if the channel count is large:
    !   S_ii | (other channels) = 1 / (S^{-1})_ii
    !   R_i^2 = 1 - S_ii|others / S_ii .
    ! Compute the ratio as (1/S_ii)/(S^{-1})_ii, avoiding product overflow.
    pure subroutine coherence_cell(s, diagonal_precision, ridge, residual, explained, r2, ok)
        complex(dp), intent(in) :: s(:,:) !! original unregularized channel CSD matrix
        complex(dp), intent(in) :: diagonal_precision(:) !! diagonal of inverted conditioned CSD
        real(dp), intent(in) :: ridge !! relative trace loading used in inversion
        real(dp), intent(out) :: residual(:) !! residual conditional spectral power
        real(dp), intent(out) :: explained(:) !! explained regularized power
        real(dp), intent(out) :: r2(:) !! squared multiple coherence
        logical, intent(out) :: ok !! all resulting values are finite and in range
        real(dp) :: mean_power, diagonal, precision, part, ratio
        integer :: n, i

        n = size(s,1)
        ok = .false.
        ! Divide before summation to reduce overflow in spectral trace.
        mean_power = 0.0_dp
        do i = 1, n
            mean_power = mean_power + real(s(i,i),dp)/real(n,dp)
        end do
        do i = 1, n
            diagonal = real(s(i,i),dp) + ridge*mean_power
            precision = real(diagonal_precision(i),dp)
            if (.not. ieee_is_finite(diagonal) .or. .not. ieee_is_finite(precision)) return
            if (diagonal <= 0.0_dp .or. precision <= 0.0_dp) return
            part = 1.0_dp/precision
            ratio = (1.0_dp/diagonal)/precision
            if (.not. ieee_is_finite(part) .or. .not. ieee_is_finite(ratio)) return
            ! Large violations indicate a non-PSD input or numerical failure.
            if (ratio > 1.0_dp + 1.0e-8_dp .or. ratio < 0.0_dp) return
            residual(i) = part
            explained(i) = max(0.0_dp,diagonal-part)
            r2(i) = max(0.0_dp,min(1.0_dp,1.0_dp-ratio))
        end do
        ok = .true.
    end subroutine coherence_cell

    function multiple_coherence_matrix(spectrum, ridge, shrinkage, auto_shrinkage, &
            effective_tapers, cv_shrinkage) result(out)
        type(multitaper_matrix_result), intent(in) :: spectrum !! Hermitian CSD (frequency,i,j)
        real(dp), intent(in), optional :: ridge !! relative diagonal loading, default 0
        real(dp), intent(in), optional :: shrinkage !! diagonal-target shrinkage in [0,1]
        logical, intent(in), optional :: auto_shrinkage !! Gaussian-Wishart plug-in shrinkage
        real(dp), intent(in), optional :: effective_tapers !! fixed-weight equivalent taper count
        logical, intent(in), optional :: cv_shrinkage !! use stored LOO shrinkage coefficients
        type(multiple_coherence_result) :: out
        type(partial_coherence_result) :: inv
        real(dp), allocatable :: residual(:), explained(:), r2(:)
        real(dp) :: lambda, nan
        integer :: nf, nc, f, i
        logical :: ok

        allocate(out%frequencies(0), out%coherence(0,0), out%conditional_power(0,0), &
            out%explained_power(0,0), out%shrinkage_used(0), out%valid(0))
        inv = partial_coherence_matrix(spectrum, ridge=ridge, shrinkage=shrinkage, &
            auto_shrinkage=auto_shrinkage, effective_tapers=effective_tapers, &
            cv_shrinkage=cv_shrinkage)
        if (inv%status == partial_coherence_invalid_input) return
        nf = size(inv%frequencies)
        nc = size(inv%precision,2)
        lambda = 0.0_dp
        if (present(ridge)) lambda = ridge
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        deallocate(out%frequencies, out%coherence, out%conditional_power, &
            out%explained_power, out%shrinkage_used, out%valid)
        allocate(out%frequencies(nf), out%coherence(nf,nc), &
            out%conditional_power(nf,nc), out%explained_power(nf,nc), &
            out%shrinkage_used(nf), out%valid(nf))
        out%frequencies = inv%frequencies
        out%shrinkage_used = inv%shrinkage_used
        out%valid = inv%valid
        out%coherence = nan
        out%conditional_power = nan
        out%explained_power = nan
        out%status = inv%status
        allocate(residual(nc), explained(nc), r2(nc))
        do f = 1, nf
            if (.not. out%valid(f)) cycle
            call coherence_cell(spectrum%power(f,:,:), &
                [(inv%precision(f,i,i),i=1,nc)],lambda,residual,explained,r2,ok)
            if (.not. ok) then
                out%valid(f) = .false.
                out%status = partial_coherence_singular
                cycle
            end if
            out%coherence(f,:) = r2
            out%conditional_power(f,:) = residual
            out%explained_power(f,:) = explained
        end do
    end function multiple_coherence_matrix

    function multiple_coherence_tf_matrix(spectrum, ridge, shrinkage, auto_shrinkage, &
            effective_tapers, cv_shrinkage) result(out)
        ! Spectral matrix cube (frequency,time,i,j).
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: spectrum
        real(dp), intent(in), optional :: ridge !! relative diagonal loading, default 0
        real(dp), intent(in), optional :: shrinkage !! diagonal-target shrinkage in [0,1]
        logical, intent(in), optional :: auto_shrinkage !! Gaussian-Wishart plug-in shrinkage
        real(dp), intent(in), optional :: effective_tapers !! fixed-weight equivalent taper count
        logical, intent(in), optional :: cv_shrinkage !! use LOO coefficients stored with spectra
        type(multiple_coherence_tf_result) :: out
        type(partial_coherence_tf_result) :: inv
        real(dp), allocatable :: residual(:), explained(:), r2(:)
        real(dp) :: lambda, nan
        integer :: nf, nt, nc, f, t, i
        logical :: ok

        allocate(out%frequencies(0), out%times(0), out%coherence(0,0,0), &
            out%conditional_power(0,0,0), out%explained_power(0,0,0), &
            out%shrinkage_used(0,0), out%valid(0,0))
        inv = partial_coherence_tf_matrix(spectrum, ridge=ridge, shrinkage=shrinkage, &
            auto_shrinkage=auto_shrinkage, effective_tapers=effective_tapers, &
            cv_shrinkage=cv_shrinkage)
        if (inv%status == partial_coherence_invalid_input) return
        nf = size(inv%frequencies)
        nt = size(inv%times)
        nc = size(inv%precision,3)
        lambda = 0.0_dp
        if (present(ridge)) lambda = ridge
        nan = ieee_value(0.0_dp,ieee_quiet_nan)
        deallocate(out%frequencies, out%times, out%coherence, out%conditional_power, &
            out%explained_power, out%shrinkage_used, out%valid)
        allocate(out%frequencies(nf), out%times(nt), out%coherence(nf,nt,nc), &
            out%conditional_power(nf,nt,nc), out%explained_power(nf,nt,nc), &
            out%shrinkage_used(nf,nt), out%valid(nf,nt))
        out%frequencies = inv%frequencies
        out%times = inv%times
        out%shrinkage_used = inv%shrinkage_used
        out%valid = inv%valid
        out%coherence = nan
        out%conditional_power = nan
        out%explained_power = nan
        out%status = inv%status
        allocate(residual(nc), explained(nc), r2(nc))
        do t = 1, nt
            do f = 1, nf
                if (.not. out%valid(f,t)) cycle
                call coherence_cell(spectrum%power(f,t,:,:), &
                    [(inv%precision(f,t,i,i),i=1,nc)],lambda,residual,explained,r2,ok)
                if (.not. ok) then
                    out%valid(f,t) = .false.
                    out%status = partial_coherence_singular
                    cycle
                end if
                out%coherence(f,t,:) = r2
                out%conditional_power(f,t,:) = residual
                out%explained_power(f,t,:) = explained
            end do
        end do
    end function multiple_coherence_tf_matrix

    function multitaper_multiple_coherence_real(x, nw, kmax, fs, nfft, detrend, &
            weighting, axis, detrend_mode, ridge, shrinkage, &
            auto_shrinkage, cv_shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! real observations, sample x channel by default
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate in Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean before tapering
        character(*), intent(in), optional :: weighting !! equal, eigen, or adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, or linear
        real(dp), intent(in), optional :: ridge !! nonnegative relative trace loading
        real(dp), intent(in), optional :: shrinkage !! explicit diagonal shrinkage [0,1]
        logical, intent(in), optional :: auto_shrinkage !! plug-in shrinkage for fixed weights
        logical, intent(in), optional :: cv_shrinkage !! LOO shrinkage including adaptive weights
        type(multiple_coherence_result) :: out
        type(multitaper_matrix_result) :: spectrum

        spectrum = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend, &
            weighting,axis,detrend_mode,cv_shrinkage)
        if (present(auto_shrinkage)) then
            if (auto_shrinkage) then
                out = multiple_coherence_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
                    auto_shrinkage=.true.,effective_tapers=spectrum%effective_tapers, &
                    cv_shrinkage=cv_shrinkage)
                return
            end if
        end if
        out = multiple_coherence_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
            cv_shrinkage=cv_shrinkage)
    end function multitaper_multiple_coherence_real

    function multitaper_multiple_coherence_complex(x, nw, kmax, fs, nfft, detrend, &
            weighting, axis, detrend_mode, ridge, shrinkage, &
            auto_shrinkage, cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex observations, sample x channel by default
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate in Hz
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean before tapering
        character(*), intent(in), optional :: weighting !! equal, eigen, or adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, or linear
        real(dp), intent(in), optional :: ridge !! nonnegative relative trace loading
        real(dp), intent(in), optional :: shrinkage !! explicit diagonal shrinkage [0,1]
        logical, intent(in), optional :: auto_shrinkage !! plug-in shrinkage for fixed weights
        logical, intent(in), optional :: cv_shrinkage !! LOO shrinkage including adaptive weights
        type(multiple_coherence_result) :: out
        type(multitaper_matrix_result) :: spectrum

        spectrum = multitaper_spectral_matrix(x,nw,kmax,fs,nfft,detrend, &
            weighting,axis,detrend_mode,cv_shrinkage)
        if (present(auto_shrinkage)) then
            if (auto_shrinkage) then
                out = multiple_coherence_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
                    auto_shrinkage=.true.,effective_tapers=spectrum%effective_tapers, &
                    cv_shrinkage=cv_shrinkage)
                return
            end if
        end if
        out = multiple_coherence_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
            cv_shrinkage=cv_shrinkage)
    end function multitaper_multiple_coherence_complex

    function multitaper_tf_multiple_coherence_real(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            auto_shrinkage,cv_shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! real observations, sample x channel by default
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! sliding-window sample count
        integer, intent(in), optional :: noverlap !! shared samples of adjacent windows
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! positive sample frequency
        integer, intent(in), optional :: nfft !! FFT length >= nperseg
        logical, intent(in), optional :: detrend !! whether to remove means
        character(*), intent(in), optional :: weighting !! equal, eigen, or adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, or linear
        real(dp), intent(in), optional :: ridge !! nonnegative relative trace loading
        real(dp), intent(in), optional :: shrinkage !! explicit diagonal shrinkage [0,1]
        logical, intent(in), optional :: auto_shrinkage !! plug-in fixed-weight shrinkage
        logical, intent(in), optional :: cv_shrinkage !! LOO shrinkage including adaptive weights
        type(multiple_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectrum

        spectrum = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,cv_shrinkage)
        if (present(auto_shrinkage)) then
            if (auto_shrinkage) then
                out = multiple_coherence_tf_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
                    auto_shrinkage=.true.,effective_tapers=spectrum%effective_tapers, &
                    cv_shrinkage=cv_shrinkage)
                return
            end if
        end if
        out = multiple_coherence_tf_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
            cv_shrinkage=cv_shrinkage)
    end function multitaper_tf_multiple_coherence_real

    function multitaper_tf_multiple_coherence_complex(x,nw,nperseg,noverlap,kmax,fs,nfft, &
            detrend,weighting,axis,detrend_mode,ridge,shrinkage, &
            auto_shrinkage,cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex observations, sample x channel by default
        real(dp), intent(in) :: nw !! positive DPSS time-bandwidth product
        integer, intent(in) :: nperseg !! sliding-window sample count
        integer, intent(in), optional :: noverlap !! shared samples of adjacent windows
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! positive sample frequency
        integer, intent(in), optional :: nfft !! FFT length >= nperseg
        logical, intent(in), optional :: detrend !! whether to remove means
        character(*), intent(in), optional :: weighting !! equal, eigen, or adaptive
        integer, intent(in), optional :: axis !! one-based sample axis, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none, constant, or linear
        real(dp), intent(in), optional :: ridge !! nonnegative relative trace loading
        real(dp), intent(in), optional :: shrinkage !! explicit diagonal shrinkage [0,1]
        logical, intent(in), optional :: auto_shrinkage !! plug-in fixed-weight shrinkage
        logical, intent(in), optional :: cv_shrinkage !! LOO shrinkage including adaptive weights
        type(multiple_coherence_tf_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: spectrum

        spectrum = multitaper_spectral_matrix_spectrogram(x,nw,nperseg,noverlap, &
            kmax,fs,nfft,detrend,weighting,axis,detrend_mode,cv_shrinkage)
        if (present(auto_shrinkage)) then
            if (auto_shrinkage) then
                out = multiple_coherence_tf_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
                    auto_shrinkage=.true.,effective_tapers=spectrum%effective_tapers, &
                    cv_shrinkage=cv_shrinkage)
                return
            end if
        end if
        out = multiple_coherence_tf_matrix(spectrum,ridge=ridge,shrinkage=shrinkage, &
            cv_shrinkage=cv_shrinkage)
    end function multitaper_tf_multiple_coherence_complex

end module scifort_signal_multiple_coherence
