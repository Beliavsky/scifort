! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Multitaper spectral estimation: SciFort extensions using Slepian tapers.
! Each orthonormal taper produces a periodogram, combined with equal,
! concentration, or Thomson-style adaptive weighting. Includes pointwise
! jackknife uncertainty for one-dimensional real and complex signals.
! Includes windowed time-frequency PSDs and the Thomson harmonic F-test.
module scifort_signal_multitaper
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, spectral_status_success, &
        spectral_status_invalid_input, fft_plan
    use scifort_signal_dpss, only : dpss, dpss_result, dpss_status_success
    implicit none
    private
    type, public :: multitaper_matrix_result
        real(dp), allocatable :: frequencies(:)
        complex(dp), allocatable :: power(:,:,:) !! frequency, channel x, channel y
        real(dp), allocatable :: cv_shrinkage(:) !! optional LOO cross-validated shrinkage (frequency)
        real(dp) :: effective_tapers = 0.0_dp !! independent fixed-weight equivalents; 0 for adaptive/invalid
        integer :: status = spectral_status_invalid_input
    end type multitaper_matrix_result

    type, public :: multitaper_channels_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: power(:,:) !! frequency, channel
        integer :: status = spectral_status_invalid_input
    end type multitaper_channels_result

    type, public :: multitaper_confidence_result
        real(dp), allocatable :: frequencies(:) !! FFT frequencies in Hz
        real(dp), allocatable :: power(:) !! PSD, using the original units
        real(dp), allocatable :: lower(:) !! log-jackknife lower confidence limit
        real(dp), allocatable :: upper(:) !! log-jackknife upper confidence limit
        real(dp), allocatable :: stderr_log(:) !! jackknife standard error of log(PSD)
        real(dp), allocatable :: equivalent_dof(:) !! approximate spectral degrees of freedom
        integer :: status = spectral_status_invalid_input
    end type multitaper_confidence_result

    ! SciFort time-frequency extension. Frequency is the first dimension;
    ! segment centre is the second; channel is the third, where present.
    type, public :: multitaper_spectrogram_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        real(dp), allocatable :: power(:,:) !! frequency, time
        integer :: status = spectral_status_invalid_input
    end type multitaper_spectrogram_result

    type, public :: multitaper_spectrogram_channels_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        real(dp), allocatable :: power(:,:,:) !! frequency, time, channel
        integer :: status = spectral_status_invalid_input
    end type multitaper_spectrogram_channels_result

    ! Cross-channel PSD convention: P_xy = conj(FFT(x))*FFT(y).
    type, public :: multitaper_spectral_matrix_spectrogram_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        complex(dp), allocatable :: power(:,:,:,:) !! frequency, time, x-channel, y-channel
        real(dp), allocatable :: cv_shrinkage(:,:) !! optional LOO shrinkage (frequency,time)
        real(dp) :: effective_tapers = 0.0_dp !! independent fixed-weight equivalents; 0 for adaptive/invalid
        integer :: status = spectral_status_invalid_input
    end type multitaper_spectral_matrix_spectrogram_result

    type, public :: multitaper_cross_spectrogram_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        complex(dp), allocatable :: power(:,:) !! P_xy(frequency,time)
        real(dp), allocatable :: auto_x(:,:) !! P_xx(frequency,time)
        real(dp), allocatable :: auto_y(:,:) !! P_yy(frequency,time)
        integer :: status = spectral_status_invalid_input
    end type multitaper_cross_spectrogram_result

    type, public :: multitaper_tf_coherence_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        real(dp), allocatable :: coherence(:,:) !! squared cross-spectrum magnitude / (P_xx*P_yy)
        integer :: status = spectral_status_invalid_input
    end type multitaper_tf_coherence_result

    type, public :: multitaper_harmonic_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: f_statistic(:)
        real(dp), allocatable :: pvalue(:)
        complex(dp), allocatable :: amplitude(:)
        integer :: status = spectral_status_invalid_input
    end type multitaper_harmonic_result

    ! Robust local background for colored-noise multitaper spectra.
    ! The ratio is descriptive; it is not a statistical p-value.
    type, public :: multitaper_noise_floor_result
        real(dp), allocatable :: frequencies(:) !! Fourier frequencies in Hz
        real(dp), allocatable :: power(:) !! original multitaper power spectrum
        real(dp), allocatable :: background(:) !! local median spectrum
        real(dp), allocatable :: whitened(:) !! power divided by local median
        integer :: status = spectral_status_invalid_input
    end type multitaper_noise_floor_result

    ! Multichannel colored-noise estimate for a single (mean) spectrum.
    type, public :: multitaper_noise_floor_channels_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: power(:,:) !! frequency, channel
        real(dp), allocatable :: background(:,:) !! frequency, channel
        real(dp), allocatable :: whitened(:,:) !! frequency, channel
        integer :: status = spectral_status_invalid_input
    end type multitaper_noise_floor_channels_result

    ! Time-frequency version of the colored-noise background estimate.
    ! Each window and channel is normalized independently along frequency.
    type, public :: multitaper_noise_floor_spectrogram_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        real(dp), allocatable :: power(:,:) !! frequency, time
        real(dp), allocatable :: background(:,:) !! frequency, time
        real(dp), allocatable :: whitened(:,:) !! frequency, time
        integer :: status = spectral_status_invalid_input
    end type multitaper_noise_floor_spectrogram_result

    type, public :: multitaper_noise_floor_spectrogram_channels_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        real(dp), allocatable :: power(:,:,:) !! frequency, time, channel
        real(dp), allocatable :: background(:,:,:) !! frequency, time, channel
        real(dp), allocatable :: whitened(:,:,:) !! frequency, time, channel
        integer :: status = spectral_status_invalid_input
    end type multitaper_noise_floor_spectrogram_channels_result

    ! Family-wise/FDR-corrected harmonic-line scan (SciFort extension).
    type, public :: multitaper_harmonic_scan_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: f_statistic(:)
        real(dp), allocatable :: pvalue(:)
        real(dp), allocatable :: adjusted_pvalue(:)
        complex(dp), allocatable :: amplitude(:)
        logical, allocatable :: significant(:)
        integer :: num_tests = 0
        integer :: status = spectral_status_invalid_input
    end type multitaper_harmonic_scan_result

    ! Incremental multitaper spectral matrices. Storage depends on window and
    ! channel count, not on the number of samples already processed.
    type, public :: multitaper_stream
        private
        integer :: nperseg = 0
        integer :: nfft = 0
        integer :: nk = 0
        integer :: hop = 0
        integer :: nc = 0
        integer :: nt = 0
        integer :: buffered = 0
        integer(int64) :: samples_seen = 0_int64
        integer(int64) :: frames_seen = 0_int64
        real(dp) :: rate = 1.0_dp
        logical :: real_input = .true.
        character(:), allocatable :: trend
        logical :: ready = .false.
        character(:), allocatable :: weighting
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: tapers(:,:), ratios(:), normalizers(:)
        real(dp), allocatable :: periodograms(:,:), weights(:,:,:), dof(:)
        complex(dp), allocatable :: window(:,:), segment(:), transform(:), spectra(:,:,:)
        complex(dp), allocatable :: summed_spectra(:,:,:)
        type(fft_plan) :: fft_state !! scratch and roots reused across windows
    contains
        procedure, public :: initialize => multitaper_stream_initialize
        procedure, public :: reset => multitaper_stream_reset
        procedure, public :: push_real => multitaper_stream_push_real
        procedure, public :: push_complex => multitaper_stream_push_complex
        procedure, public :: push_real_channels => multitaper_stream_push_real_channels
        procedure, public :: push_complex_channels => multitaper_stream_push_complex_channels
        generic, public :: push => push_real, push_complex, push_real_channels, push_complex_channels
        procedure, public :: push_noise_floor_real => multitaper_stream_floor_real
        procedure, public :: push_noise_floor_complex => multitaper_stream_floor_complex
        procedure, public :: push_noise_floor_real_channels => multitaper_stream_floor_real_channels
        procedure, public :: push_noise_floor_complex_channels => multitaper_stream_floor_complex_channels
        generic, public :: push_noise_floor => push_noise_floor_real, push_noise_floor_complex, &
            push_noise_floor_real_channels, push_noise_floor_complex_channels
        procedure, public :: mean_noise_floor => multitaper_stream_mean_noise_floor
        procedure, public :: mean_matrix => multitaper_stream_mean_matrix
        procedure, public :: sample_count => multitaper_stream_sample_count
        procedure, public :: window_count => multitaper_stream_window_count
        procedure, private :: output_shell => multitaper_stream_output_shell
        procedure, private :: emit => multitaper_stream_emit
        procedure, private :: compute_window => multitaper_stream_compute_window
    end type multitaper_stream

    interface multitaper_spectrogram
        module procedure multitaper_spectrogram_real
        module procedure multitaper_spectrogram_complex
        module procedure multitaper_spectrogram_real_channels
        module procedure multitaper_spectrogram_complex_channels
    end interface
    interface multitaper_spectral_matrix_spectrogram
        module procedure multitaper_matrix_spectrogram_real
        module procedure multitaper_matrix_spectrogram_complex
    end interface
    interface multitaper_cross_spectrogram
        module procedure multitaper_cross_spectrogram_real
        module procedure multitaper_cross_spectrogram_complex
        module procedure multitaper_cross_spectrogram_real_complex
        module procedure multitaper_cross_spectrogram_complex_real
    end interface
    interface multitaper_tf_coherence
        module procedure multitaper_tf_coherence_real
        module procedure multitaper_tf_coherence_complex
        module procedure multitaper_tf_coherence_real_complex
        module procedure multitaper_tf_coherence_complex_real
    end interface
    interface multitaper_harmonic_ftest
        module procedure multitaper_harmonic_ftest_real
        module procedure multitaper_harmonic_ftest_complex
    end interface
    interface multitaper_harmonic_scan
        module procedure multitaper_harmonic_scan_real
        module procedure multitaper_harmonic_scan_complex
    end interface
    interface multitaper_noise_floor
        module procedure multitaper_noise_floor_real
        module procedure multitaper_noise_floor_complex
    end interface
    interface multitaper_noise_floor_spectrogram
        module procedure multitaper_noise_floor_spectrogram_real
        module procedure multitaper_noise_floor_spectrogram_complex
        module procedure multitaper_noise_floor_spectrogram_real_channels
        module procedure multitaper_noise_floor_spectrogram_complex_channels
    end interface
    interface multitaper_psd
        module procedure multitaper_psd_real
        module procedure multitaper_psd_complex
        module procedure multitaper_psd_real_channels
        module procedure multitaper_psd_complex_channels
    end interface
    interface multitaper_spectral_matrix
        module procedure multitaper_matrix_real
        module procedure multitaper_matrix_complex
    end interface
    interface multitaper_csd
        module procedure multitaper_csd_real
        module procedure multitaper_csd_complex
        module procedure multitaper_csd_real_complex
        module procedure multitaper_csd_complex_real
    end interface
    interface multitaper_coherence
        module procedure multitaper_coherence_real
        module procedure multitaper_coherence_complex
    end interface
    interface multitaper_psd_ci
        module procedure multitaper_psd_ci_real
        module procedure multitaper_psd_ci_complex
    end interface
    public :: multitaper_psd, multitaper_spectral_matrix, multitaper_csd, multitaper_coherence
    public :: multitaper_effective_tapers
    public :: multitaper_psd_ci, multitaper_spectrogram, multitaper_harmonic_ftest
    public :: multitaper_harmonic_scan
    public :: multitaper_noise_floor
    public :: multitaper_noise_floor_spectrogram
    public :: multitaper_spectral_matrix_spectrogram, multitaper_cross_spectrogram, multitaper_tf_coherence
contains
    ! Remove a per-segment least-squares polynomial in centered coordinates.
    ! Removing the mean before the slope avoids cancellation from a large
    ! offset. The same complex slope calculation handles real input too.
    subroutine detrend_multitaper_segment(segment,trend)
        complex(dp), intent(inout) :: segment(:) !! samples in a time window
        character(*), intent(in) :: trend !! none, constant or linear
        real(dp) :: coordinate, denominator, midpoint
        complex(dp) :: slope
        integer :: i, n

        if (trend == 'none') return
        n = size(segment)
        segment = segment - sum(segment)/real(n,dp)
        if (trend /= 'linear' .or. n <= 1) return
        midpoint = 0.5_dp*real(n+1,dp)
        denominator = real(n,dp)*(real(n,dp)**2-1.0_dp)/12.0_dp
        slope = (0.0_dp,0.0_dp)
        do i = 1,n
            coordinate = real(i,dp) - midpoint
            slope = slope + coordinate*segment(i)
        end do
        slope = slope/denominator
        do i = 1,n
            coordinate = real(i,dp) - midpoint
            segment(i) = segment(i) - coordinate*slope
        end do
    end subroutine detrend_multitaper_segment

    ! Frequency-dependent Thomson adaptive taper weights. The input periodograms
    ! must use two-sided density units (before doubling real non-DC bins).
    ! Returned q(t,f) is nonnegative with sum(q(:,f)) = 1, so the routine
    ! also works for equal/eigen weights and provides an effective DOF.
    ! Fixed deterministic taper weights in a weighted complex-Wishart model.
    ! Adaptive weights depend on data and are explicitly not supported.
    pure function multitaper_effective_tapers(weights) result(n)
        real(dp), intent(in) :: weights(:) !! finite nonnegative fixed weights, not all zero
        real(dp) :: n, scale, total
        real(dp), allocatable :: normalized(:)

        n = 0.0_dp
        if (size(weights) < 1) return
        if (.not. all(ieee_is_finite(weights))) return
        if (any(weights < 0.0_dp)) return
        scale = maxval(weights)
        if (scale <= 0.0_dp) return
        normalized = weights/scale
        total = sum(normalized)
        n = (total/sqrt(sum(normalized**2)))**2
    end function multitaper_effective_tapers

    subroutine multitaper_weights(periodograms, eigenvalues, variance, method, q, dof)
        real(dp), intent(in) :: periodograms(:,:) !! PSD by taper, frequency
        real(dp), intent(in) :: eigenvalues(:) !! DPSS spectral concentration
        real(dp), intent(in) :: variance !! two-sided noise-floor spectral density
        character(*), intent(in) :: method !! equal, eigen or adaptive
        real(dp), intent(out) :: q(:,:) !! normalized squared weights
        real(dp), intent(out) :: dof(:) !! 2/sum(q**2)
        real(dp) :: s, snext, w, denominator, total, cutoff, r, lambda
        integer :: nf, nt, f, t, iter

        nt = size(periodograms,1)
        nf = size(periodograms,2)
        q = 0.0_dp
        do f = 1, nf
            if (method == 'equal') then
                q(:,f) = 1.0_dp/real(nt,dp)
            else if (method == 'eigen' .or. nt < 3) then
                q(:,f) = max(eigenvalues,0.0_dp)
                total = sum(q(:,f))
                if (total > 0.0_dp) then
                    q(:,f) = q(:,f)/total
                else
                    q(:,f) = 1.0_dp/real(nt,dp)
                end if
            else
                ! Start from the first two (best-concentrated) tapers.
                total = sum(max(eigenvalues(1:2),0.0_dp))
                s = 0.0_dp
                if (total > 0.0_dp) then
                    s = sum(max(eigenvalues(1:2),0.0_dp)*periodograms(1:2,f))/total
                end if
                cutoff = max(tiny(1.0_dp)*1000.0_dp, epsilon(1.0_dp)*max(s,variance))
                do iter = 1, 200
                    total = 0.0_dp
                    do t = 1, nt
                        lambda = max(0.0_dp,min(1.0_dp,eigenvalues(t)))
                        denominator = lambda*s+(1.0_dp-lambda)*variance
                        w = 0.0_dp
                        if (denominator > cutoff) w = lambda*(s/denominator)**2
                        q(t,f) = w
                        total = total + w
                    end do
                    if (total <= tiny(1.0_dp)) then
                        q(:,f) = max(eigenvalues,0.0_dp)
                        total = sum(q(:,f))
                        if (total <= tiny(1.0_dp)) then
                            q(:,f) = 1.0_dp/real(nt,dp)
                        else
                            q(:,f) = q(:,f)/total
                        end if
                        exit
                    end if
                    q(:,f) = q(:,f)/total
                    snext = sum(q(:,f)*periodograms(:,f))
                    r = abs(snext-s)/max(max(snext,s),cutoff)
                    s = snext
                    if (r < 1.0e-10_dp) exit
                end do
            end if
            dof(f) = 2.0_dp/max(sum(q(:,f)**2),tiny(1.0_dp))
        end do
    end subroutine multitaper_weights

    ! Cross-validated diagonal-target shrinkage from individual tapered FFTs.
    ! For each omitted taper k, reconstruct S_{-k} using weights refitted
    ! without its periodogram; compare its off-diagonals with the held-out
    ! outer product H_k. The risk-minimizing real slope for this prediction
    ! loss is beta=sum Re(conjg(S_{-k,ij})*H_{k,ij}) /
    !             sum |S_{-k,ij}|**2, and alpha=clip(1-beta,0,1).
    ! For adaptive weighting, the variance used for Thomson weights is the
    ! original segment variance, but the frequency-local weights are refit
    ! without the held-out taper. This is a predictive CV criterion, not an
    ! unbiased estimator of the risk of the full weighted cross spectrum.
    subroutine multitaper_cv_diagonal(spectra, eigenvalues, variances, method, alpha)
        complex(dp), intent(in) :: spectra(:,:,:) !! normalized complex FFTs (frequency,taper,channel)
        real(dp), intent(in) :: eigenvalues(:) !! DPSS concentration, one per taper
        real(dp), intent(in) :: variances(:) !! two-sided segment variances by channel
        character(*), intent(in) :: method !! equal, eigen or adaptive
        real(dp), intent(out) :: alpha(:) !! [0,1] diagonal shrinkage by frequency

        real(dp), allocatable :: periods(:,:), weights(:,:,:), dof(:), numerator(:), denominator(:)
        complex(dp) :: train, validation
        real(dp) :: magnitude
        integer :: nf, nt, nc, omitted, c, i, j, f, t, r

        nf = size(spectra,1)
        nt = size(spectra,2)
        nc = size(spectra,3)
        allocate(periods(nt-1,nf),weights(nt-1,nf,nc),dof(nf))
        allocate(numerator(nf),denominator(nf))
        numerator = 0.0_dp
        denominator = 0.0_dp
        do omitted = 1, nt
            do c = 1, nc
                r = 0
                do t = 1, nt
                    if (t == omitted) cycle
                    r = r+1
                    periods(r,:) = abs(spectra(:,t,c))**2
                end do
                call multitaper_weights(periods, &
                    [eigenvalues(:omitted-1),eigenvalues(omitted+1:)], &
                    variances(c),method,weights(:,:,c),dof)
            end do
            do f = 1, nf
                magnitude = maxval(abs(spectra(f,:,:)))
                if (magnitude <= 0.0_dp) cycle
                do j = 2, nc
                    do i = 1, j-1
                        train = (0.0_dp,0.0_dp)
                        r = 0
                        do t = 1, nt
                            if (t == omitted) cycle
                            r = r+1
                            train = train + sqrt(weights(r,f,i)*weights(r,f,j))* &
                                conjg(spectra(f,t,i)/magnitude)*(spectra(f,t,j)/magnitude)
                        end do
                        validation = conjg(spectra(f,omitted,i)/magnitude)* &
                            (spectra(f,omitted,j)/magnitude)
                        numerator(f) = numerator(f)+real(conjg(train)*validation,dp)
                        denominator(f) = denominator(f)+abs(train)**2
                    end do
                end do
            end do
        end do
        ! When every training off-diagonal vanishes, shrink fully. At an
        ! all-zero bin the subsequent Cholesky inversion marks it singular.
        alpha = 1.0_dp
        where (denominator > tiny(1.0_dp))
            alpha = max(0.0_dp,min(1.0_dp,1.0_dp-numerator/max(denominator,tiny(1.0_dp))))
        end where
    end subroutine multitaper_cv_diagonal


    ! The regularized incomplete beta is used only to compute confidence
    ! multipliers.  The continued fraction avoids tables and external deps.
    pure real(dp) function beta_contfrac(a,b,x) result(h)
        real(dp), intent(in) :: a,b,x
        real(dp) :: aa, c, d, del, qab, qam, qap, pmin
        integer :: m, m2
        pmin = tiny(1.0_dp)/epsilon(1.0_dp)
        qab = a+b
        qap = a+1.0_dp
        qam = a-1.0_dp
        c = 1.0_dp
        d = 1.0_dp-qab*x/qap
        if (abs(d) < pmin) d = pmin
        d = 1.0_dp/d
        h = d
        do m = 1, 256
            m2 = 2*m
            aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
            d = 1.0_dp+aa*d
            if (abs(d) < pmin) d = pmin
            c = 1.0_dp+aa/c
            if (abs(c) < pmin) c = pmin
            d = 1.0_dp/d
            h = h*d*c
            aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
                ((a+real(m2,dp))*(qap+real(m2,dp)))
            d = 1.0_dp+aa*d
            if (abs(d) < pmin) d = pmin
            c = 1.0_dp+aa/c
            if (abs(c) < pmin) c = pmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= 8.0_dp*epsilon(1.0_dp)) exit
        end do
    end function beta_contfrac

    pure real(dp) function regularized_beta(a,b,x) result(value)
        real(dp), intent(in) :: a,b,x
        real(dp) :: bt
        if (x <= 0.0_dp) then
            value = 0.0_dp
        else if (x >= 1.0_dp) then
            value = 1.0_dp
        else
            bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+ &
                a*log(x)+b*log(1.0_dp-x))
            if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
                value = bt*beta_contfrac(a,b,x)/a
            else
                value = 1.0_dp-bt*beta_contfrac(b,a,1.0_dp-x)/b
            end if
            value = max(0.0_dp,min(1.0_dp,value))
        end if
    end function regularized_beta

    pure real(dp) function student_t_critical(df, confidence) result(critical)
        integer, intent(in) :: df !! K-1 jackknife degrees of freedom
        real(dp), intent(in) :: confidence !! two-sided coverage probability
        real(dp) :: hi, lo, mid, tail, nu
        integer :: i
        nu = real(df,dp)
        lo = 0.0_dp
        hi = 1.0_dp
        do i = 1, 100
            tail = 0.5_dp*regularized_beta(0.5_dp*nu,0.5_dp,nu/(nu+hi*hi))
            if (tail <= 0.5_dp*(1.0_dp-confidence)) exit
            hi = 2.0_dp*hi
        end do
        do i = 1, 65
            mid = 0.5_dp*(lo+hi)
            tail = 0.5_dp*regularized_beta(0.5_dp*nu,0.5_dp,nu/(nu+mid*mid))
            if (tail > 0.5_dp*(1.0_dp-confidence)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        critical = 0.5_dp*(lo+hi)
    end function student_t_critical

    function multitaper_ci_core(x, real_input, nw, kmax, fs, nfft, detrend, weighting, confidence,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! signal
        logical, intent(in) :: real_input !! one-sided spectrum for real data
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! taper count, must be >= 2
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! subtract mean
        character(*), intent(in), optional :: weighting !! equal, eigen or adaptive
        real(dp), intent(in), optional :: confidence !! coverage, default 0.95
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_confidence_result) :: out
        type(dpss_result) :: windows
        type(fft_plan) :: fft_state
        integer :: fft_status
        complex(dp), allocatable :: signal(:), fft(:)
        real(dp), allocatable :: raw(:,:), weights(:,:), dof(:), partial(:,:), rweights(:,:), rdof(:), jklog(:)
        real(dp) :: rate, coverage, variance, factor, crit, ps, theta, v, lowlog, uplog
        real(dp) :: s, meanlog, logmin, logmax
        integer :: n, nf, nk, nt, j, f, omit, ibin, sb
        character(:), allocatable :: combine, trend

        allocate(out%frequencies(0),out%power(0),out%lower(0),out%upper(0), &
            out%stderr_log(0),out%equivalent_dof(0))
        n = size(x)
        if (n < 2 .or. .not. ieee_is_finite(nw)) return
        if (.not. all(ieee_is_finite(real(x,dp))) .or. .not. all(ieee_is_finite(aimag(x)))) return
        if (nw <= 0.0_dp .or. nw >= 0.5_dp*real(n,dp)) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate) .or. rate <= 0.0_dp) return
        nf = n
        if (present(nfft)) nf = nfft
        if (nf < n) return
        nt = min(n,max(1,int(2.0_dp*nw)-1))
        if (present(kmax)) nt = kmax
        if (nt < 2 .or. nt > n) return
        coverage = 0.95_dp
        if (present(confidence)) coverage = confidence
        if (.not. ieee_is_finite(coverage)) return
        if (coverage <= 0.0_dp .or. coverage >= 1.0_dp) return
        combine = 'equal'
        if (present(weighting)) combine = trim(weighting)
        if (combine /= 'equal' .and. combine /= 'eigen' .and. combine /= 'adaptive') return
        trend = 'constant'
        if (present(detrend)) then
            if (detrend) then
                trend = 'constant'
            else
                trend = 'none'
            end if
        end if
        if (present(detrend_mode)) trend = trim(detrend_mode)
        if (trend /= 'none' .and. trend /= 'constant' .and. trend /= 'linear') return
        windows = dpss(n,nw,kmax=nt,sym=.false.,norm='l2')
        if (windows%status /= dpss_status_success) return
        call fft_state%initialize(nf,fft_status)
        if (fft_status /= spectral_status_success) return
        nk = nf
        if (real_input) nk = nf/2+1
        allocate(signal(n),fft(nf),raw(nt,nk),weights(nt,nk),dof(nk), &
            partial(nt-1,nk),rweights(nt-1,nk),rdof(nk),jklog(nt))
        signal = x
        call detrend_multitaper_segment(signal,trend)
        variance = sum(abs(signal)**2)/real(n,dp)/rate
        do j = 1, nt
            fft = (0.0_dp,0.0_dp)
            fft(1:n) = signal*windows%tapers(:,j)
            call fft_state%execute(fft,fft_status)
            factor = rate*sum(windows%tapers(:,j)**2)
            raw(j,:) = abs(fft(1:nk))**2/factor
        end do
        call multitaper_weights(raw,windows%ratios,variance,combine,weights,dof)
        deallocate(out%frequencies,out%power,out%lower,out%upper,out%stderr_log,out%equivalent_dof)
        allocate(out%frequencies(nk),out%power(nk),out%lower(nk),out%upper(nk), &
            out%stderr_log(nk),out%equivalent_dof(nk))
        out%equivalent_dof = dof
        crit = student_t_critical(nt-1,coverage)
        logmin = log(tiny(1.0_dp))
        logmax = log(huge(1.0_dp))-1.0_dp
        do f = 1, nk
            ibin = f-1
            sb = ibin
            if (.not. real_input .and. ibin > (nf-1)/2) sb = ibin-nf
            out%frequencies(f) = real(sb,dp)*rate/real(nf,dp)
            factor = 1.0_dp
            if (real_input .and. f > 1) then
                factor = 2.0_dp
                if (mod(nf,2) == 0 .and. f == nk) factor = 1.0_dp
            end if
            ps = sum(weights(:,f)*raw(:,f))
            out%power(f) = factor*ps
            if (ps == 0.0_dp) then
                out%stderr_log(f) = 0.0_dp
                out%lower(f) = 0.0_dp
                out%upper(f) = 0.0_dp
                cycle
            end if
            do omit = 1, nt
                if (omit > 1) partial(1:omit-1,f) = raw(1:omit-1,f)
                if (omit < nt) partial(omit:nt-1,f) = raw(omit+1:nt,f)
                ! Refit adaptive weights after deleting each taper; for
                ! equal/eigen, recomputing simply renormalizes weights.
                if (omit > 1 .and. omit < nt) then
                    call multitaper_weights(partial(1:nt-1,f:f), &
                        [windows%ratios(1:omit-1),windows%ratios(omit+1:nt)], &
                        variance,combine,rweights(:,f:f),rdof(f:f))
                else if (omit == 1) then
                    call multitaper_weights(partial(1:nt-1,f:f),windows%ratios(2:nt), &
                        variance,combine,rweights(:,f:f),rdof(f:f))
                else
                    call multitaper_weights(partial(1:nt-1,f:f),windows%ratios(1:nt-1), &
                        variance,combine,rweights(:,f:f),rdof(f:f))
                end if
                s = sum(rweights(:,f)*partial(1:nt-1,f))
                jklog(omit) = log(max(s,tiny(1.0_dp)))
            end do
            meanlog = sum(jklog)/real(nt,dp)
            v = real(nt-1,dp)/real(nt,dp)*sum((jklog-meanlog)**2)
            out%stderr_log(f) = sqrt(max(0.0_dp,v))
            theta = log(ps)
            lowlog = max(logmin,theta-crit*out%stderr_log(f))
            uplog = min(logmax,theta+crit*out%stderr_log(f))
            out%lower(f) = factor*exp(lowlog)
            out%upper(f) = factor*exp(uplog)
        end do
        out%status = spectral_status_success
    end function multitaper_ci_core

    function multitaper_psd_ci_real(x,nw,kmax,fs,nfft,detrend,weighting,confidence,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! real-valued samples
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! tapers (at least 2)
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! subtract mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        real(dp), intent(in), optional :: confidence !! coverage, default 0.95
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_confidence_result) :: out
        out = multitaper_ci_core(cmplx(x,0.0_dp,dp),.true.,nw,kmax,fs,nfft,detrend,weighting,confidence,detrend_mode)
    end function multitaper_psd_ci_real

    function multitaper_psd_ci_complex(x,nw,kmax,fs,nfft,detrend,weighting,confidence,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! complex-valued samples
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! tapers (at least 2)
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! subtract mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        real(dp), intent(in), optional :: confidence !! coverage, default 0.95
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_confidence_result) :: out
        out = multitaper_ci_core(x,.false.,nw,kmax,fs,nfft,detrend,weighting,confidence,detrend_mode)
    end function multitaper_psd_ci_complex

    function multitaper_psd_real(x, nw, kmax, fs, nfft, detrend, weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! real signal, one-dimensional
        real(dp), intent(in) :: nw !! time-half-bandwidth product (0 < NW < N/2)
        integer, intent(in), optional :: kmax !! default max(1, int(2*NW)-1)
        real(dp), intent(in), optional :: fs !! sampling frequency, default 1
        integer, intent(in), optional :: nfft !! number of FFT points, default N
        logical, intent(in), optional :: detrend !! subtract arithmetic mean, default true
        character(*), intent(in), optional :: weighting !! 'equal' (default) or 'eigen'
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(spectral_result) :: out

        type(dpss_result) :: tapers
        type(fft_plan) :: fft_state
        integer :: fft_status
        complex(dp), allocatable :: spectrum(:)
        real(dp), allocatable :: xcenter(:)
        complex(dp), allocatable :: detrended(:)
        real(dp) :: rate, wt, weightsum, factor
        integer :: n, nf, k, j, nk, i
        character(:), allocatable :: combine, trend
        type(multitaper_channels_result) :: channels

        allocate(out%frequencies(0),out%power(0))
        n = size(x)
        if (n < 2 .or. .not. ieee_is_finite(nw)) return
        if (any(.not. ieee_is_finite(x))) return
        if (nw <= 0.0_dp .or. nw >= 0.5_dp*real(n,dp)) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp) return
        nf = n
        if (present(nfft)) nf = nfft
        if (nf < n) return
        k = max(1, int(2.0_dp*nw)-1)
        k = min(k,n)
        if (present(kmax)) k = kmax
        if (k < 1 .or. k > n) return
        combine = 'equal'
        if (present(weighting)) combine = trim(weighting)
        if (combine /= 'equal' .and. combine /= 'eigen' .and. combine /= 'adaptive') return
        if (combine == 'adaptive') then
            channels = multitaper_channels_core(cmplx(reshape(x,[n,1]),0.0_dp,dp), &
                .true.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
            out%frequencies = channels%frequencies
            if (channels%status == spectral_status_success) then
                out%power = channels%power(:,1)
            end if
            out%status = channels%status
            return
        end if
        trend = 'constant'
        if (present(detrend)) then
            if (detrend) then
                trend = 'constant'
            else
                trend = 'none'
            end if
        end if
        if (present(detrend_mode)) trend = trim(detrend_mode)
        if (trend /= 'none' .and. trend /= 'constant' .and. trend /= 'linear') return
        tapers = dpss(n,nw,kmax=k,sym=.false.,norm='l2')
        if (tapers%status /= dpss_status_success) return
        call fft_state%initialize(nf,fft_status)
        if (fft_status /= spectral_status_success) return

        deallocate(out%frequencies,out%power)
        nk = nf/2 + 1
        allocate(out%frequencies(nk),out%power(nk))
        allocate(xcenter(n),spectrum(nf))
        out%power = 0.0_dp
        do i = 1, nk
            out%frequencies(i) = real(i-1,dp)*rate/real(nf,dp)
        end do
        xcenter = x
        if (trend /= 'none') then
            detrended = cmplx(x,0.0_dp,dp)
            call detrend_multitaper_segment(detrended,trend)
            xcenter = real(detrended,dp)
        end if
        weightsum = 0.0_dp
        do j = 1, k
            wt = 1.0_dp
            if (combine == 'eigen') wt = tapers%ratios(j)
            if (wt <= 0.0_dp) cycle
            spectrum = (0.0_dp,0.0_dp)
            spectrum(1:n) = cmplx(xcenter*tapers%tapers(:,j),0.0_dp,dp)
            call fft_state%execute(spectrum,fft_status)
            factor = rate * sum(tapers%tapers(:,j)**2)
            out%power = out%power + wt * abs(spectrum(1:nk))**2 / factor
            weightsum = weightsum + wt
        end do
        if (weightsum <= 0.0_dp) then
            out%status = spectral_status_invalid_input
            return
        end if
        out%power = out%power / weightsum
        if (mod(nf,2) == 0) then
            if (nk > 2) out%power(2:nk-1) = 2.0_dp*out%power(2:nk-1)
        else
            if (nk > 1) out%power(2:nk) = 2.0_dp*out%power(2:nk)
        end if
        out%status = spectral_status_success
    end function multitaper_psd_real

    ! Shared spectral matrix kernel. Samples are rows and channels are columns.
    ! Conjugating the first FFT follows scipy.signal.csd's Pxy convention.
    function multitaper_matrix_core(x, real_input, nw, kmax, fs, nfft, detrend, weighting,detrend_mode,cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! samples by channels
        logical, intent(in) :: real_input !! return one-sided for purely real data
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! tapers to average
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length, at least number of samples
        logical, intent(in), optional :: detrend !! remove each channel mean (default true)
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        logical, intent(in), optional :: cv_shrinkage !! leave-one-taper-out shrinkage coefficients
        type(multitaper_matrix_result) :: out

        type(dpss_result) :: windows
        type(fft_plan) :: fft_state
        integer :: fft_status
        complex(dp), allocatable :: samples(:,:), transforms(:,:), allfft(:,:,:)
        real(dp), allocatable :: pw(:,:), q(:,:,:), dof(:), variance(:)
        real(dp) :: rate, wt, weight_sum, factor
        integer :: n, nc, nf, nk, nt, t, c, d, bin, signed_bin
        logical :: cross_validate
        character(:), allocatable :: combine, trend

        allocate(out%frequencies(0),out%power(0,0,0))
        n = size(x,1)
        nc = size(x,2)
        if (n < 2 .or. nc < 1) return
        if (.not. ieee_is_finite(nw)) return
        if (nw <= 0.0_dp .or. nw >= real(n,dp)/2.0_dp) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp) return
        nf = n
        if (present(nfft)) nf = nfft
        if (nf < n) return
        nt = min(n,max(1,int(2.0_dp*nw)-1))
        if (present(kmax)) nt = kmax
        if (nt < 1 .or. nt > n) return
        combine = 'equal'
        if (present(weighting)) combine = trim(weighting)
        if (combine /= 'equal' .and. combine /= 'eigen' .and. combine /= 'adaptive') return
        cross_validate = .false.
        if (present(cv_shrinkage)) cross_validate = cv_shrinkage
        if (cross_validate .and. (nt < 3 .or. nc < 2)) return
        trend = 'constant'
        if (present(detrend)) then
            if (detrend) then
                trend = 'constant'
            else
                trend = 'none'
            end if
        end if
        if (present(detrend_mode)) trend = trim(detrend_mode)
        if (trend /= 'none' .and. trend /= 'constant' .and. trend /= 'linear') return
        windows = dpss(n,nw,kmax=nt,sym=.false.,norm='l2')
        if (windows%status /= dpss_status_success) return
        call fft_state%initialize(nf,fft_status)
        if (fft_status /= spectral_status_success) return

        nk = nf
        if (real_input) nk = nf/2 + 1
        deallocate(out%frequencies,out%power)
        allocate(out%frequencies(nk),out%power(nk,nc,nc))
        allocate(samples(n,nc),transforms(nf,nc))
        samples = x
        do c = 1, nc
            call detrend_multitaper_segment(samples(:,c),trend)
        end do
        do bin = 0, nk-1
            signed_bin = bin
            if (.not. real_input .and. bin > (nf-1)/2) signed_bin = bin-nf
            out%frequencies(bin+1) = real(signed_bin,dp)*rate/real(nf,dp)
        end do
        out%power = (0.0_dp,0.0_dp)
        if (combine == 'adaptive' .or. cross_validate) then
            allocate(allfft(nk,nt,nc),pw(nt,nk),q(nt,nk,nc),dof(nk),variance(nc))
            do c = 1, nc
                variance(c) = sum(abs(samples(:,c))**2)/real(n,dp)/rate
                do t = 1, nt
                    transforms(:,c) = (0.0_dp,0.0_dp)
                    transforms(1:n,c) = samples(:,c)*windows%tapers(:,t)
                    call fft_state%execute(transforms(:,c),fft_status)
                    allfft(:,t,c) = transforms(1:nk,c) / &
                        sqrt(rate*sum(windows%tapers(:,t)**2))
                    pw(t,:) = abs(allfft(:,t,c))**2
                end do
                call multitaper_weights(pw,windows%ratios,variance(c),combine,q(:,:,c),dof)
            end do
            do d = 1, nc
                do c = 1, nc
                    do t = 1, nt
                        out%power(:,c,d) = out%power(:,c,d) + &
                            sqrt(q(t,:,c)*q(t,:,d))*conjg(allfft(:,t,c))*allfft(:,t,d)
                    end do
                end do
            end do
            if (real_input) then
                if (mod(nf,2) == 0) then
                    if (nk > 2) out%power(2:nk-1,:,:) = 2.0_dp*out%power(2:nk-1,:,:)
                else
                    if (nk > 1) out%power(2:nk,:,:) = 2.0_dp*out%power(2:nk,:,:)
                end if
            end if
            if (cross_validate) then
                allocate(out%cv_shrinkage(nk))
                call multitaper_cv_diagonal(allfft,windows%ratios,variance,combine,out%cv_shrinkage)
            end if
            if (combine == 'equal') out%effective_tapers = real(nt,dp)
            if (combine == 'eigen') out%effective_tapers = &
                multitaper_effective_tapers(max(windows%ratios,0.0_dp))
            out%status = spectral_status_success
            return
        end if
        weight_sum = 0.0_dp
        do t = 1, nt
            wt = 1.0_dp
            if (combine == 'eigen') wt = windows%ratios(t)
            if (wt <= 0.0_dp) cycle
            transforms = (0.0_dp,0.0_dp)
            do c = 1, nc
                transforms(1:n,c) = samples(:,c)*windows%tapers(:,t)
                call fft_state%execute(transforms(:,c),fft_status)
            end do
            factor = wt/(rate*sum(windows%tapers(:,t)**2))
            do d = 1, nc
                do c = 1, nc
                    out%power(:,c,d) = out%power(:,c,d) + &
                        factor*conjg(transforms(1:nk,c))*transforms(1:nk,d)
                end do
            end do
            weight_sum = weight_sum + wt
        end do
        if (weight_sum <= 0.0_dp) then
            out%status = spectral_status_invalid_input
            return
        end if
        out%power = out%power/weight_sum
        if (real_input) then
            if (mod(nf,2) == 0) then
                if (nk > 2) out%power(2:nk-1,:,:) = 2.0_dp*out%power(2:nk-1,:,:)
            else
                if (nk > 1) out%power(2:nk,:,:) = 2.0_dp*out%power(2:nk,:,:)
            end if
        end if
        if (combine == 'equal') out%effective_tapers = real(nt,dp)
        if (combine == 'eigen') out%effective_tapers = &
            multitaper_effective_tapers(max(windows%ratios,0.0_dp))
        out%status = spectral_status_success
    end function multitaper_matrix_core

    function multitaper_matrix_real(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,cv_shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! real input matrix
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove channel means
        character(*), intent(in), optional :: weighting !! equal or eigen
        integer, intent(in), optional :: axis !! time axis: 1 (default) or 2
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        logical, intent(in), optional :: cv_shrinkage !! compute LOO coefficients when true
        type(multitaper_matrix_result) :: out
        integer :: ax
        ax = 1
        if (present(axis)) ax = axis
        select case (ax)
        case (1)
            out = multitaper_matrix_core(cmplx(x,0.0_dp,dp),.true.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode,cv_shrinkage)
        case (2)
            out = multitaper_matrix_core(cmplx(transpose(x),0.0_dp,dp), &
                .true.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode,cv_shrinkage)
        case default
            allocate(out%frequencies(0),out%power(0,0,0))
        end select
    end function multitaper_matrix_real

    function multitaper_matrix_complex(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode,cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex input matrix
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove channel means
        character(*), intent(in), optional :: weighting !! equal or eigen
        integer, intent(in), optional :: axis !! time axis: 1 (default) or 2
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        logical, intent(in), optional :: cv_shrinkage !! compute LOO coefficients when true
        type(multitaper_matrix_result) :: out
        integer :: ax
        ax = 1
        if (present(axis)) ax = axis
        select case (ax)
        case (1)
            out = multitaper_matrix_core(x,.false.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode,cv_shrinkage)
        case (2)
            out = multitaper_matrix_core(transpose(x),.false.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode,cv_shrinkage)
        case default
            allocate(out%frequencies(0),out%power(0,0,0))
        end select
    end function multitaper_matrix_complex


    ! Compute the diagonal spectra directly, without allocating an
    ! O(nfrequency * nchannels**2) cross-spectral matrix.
    function multitaper_channels_core(x,real_input,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:,:) !! samples by channels
        logical, intent(in) :: real_input !! use one-sided result for real input
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove per-channel mean
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_channels_result) :: out

        type(dpss_result) :: windows
        type(fft_plan) :: fft_state
        integer :: fft_status
        complex(dp), allocatable :: samples(:,:), transforms(:)
        real(dp), allocatable :: pw(:,:), q(:,:), dof(:)
        real(dp) :: rate, wt, weight_sum, factor, variance
        integer :: n, nc, nf, nk, nt, t, c, bin, signed_bin
        character(:), allocatable :: combine, trend

        allocate(out%frequencies(0),out%power(0,0))
        n = size(x,1)
        nc = size(x,2)
        if (n < 2 .or. nc < 1) return
        if (.not. ieee_is_finite(nw)) return
        if (nw <= 0.0_dp .or. nw >= real(n,dp)/2.0_dp) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp) return
        nf = n
        if (present(nfft)) nf = nfft
        if (nf < n) return
        nt = min(n,max(1,int(2.0_dp*nw)-1))
        if (present(kmax)) nt = kmax
        if (nt < 1 .or. nt > n) return
        combine = 'equal'
        if (present(weighting)) combine = trim(weighting)
        if (combine /= 'equal' .and. combine /= 'eigen' .and. combine /= 'adaptive') return
        trend = 'constant'
        if (present(detrend)) then
            if (detrend) then
                trend = 'constant'
            else
                trend = 'none'
            end if
        end if
        if (present(detrend_mode)) trend = trim(detrend_mode)
        if (trend /= 'none' .and. trend /= 'constant' .and. trend /= 'linear') return
        windows = dpss(n,nw,kmax=nt,sym=.false.,norm='l2')
        if (windows%status /= dpss_status_success) return
        call fft_state%initialize(nf,fft_status)
        if (fft_status /= spectral_status_success) return

        nk = nf
        if (real_input) nk = nf/2+1
        deallocate(out%frequencies,out%power)
        allocate(out%frequencies(nk),out%power(nk,nc))
        allocate(samples(n,nc),transforms(nf))
        samples = x
        do c = 1, nc
            call detrend_multitaper_segment(samples(:,c),trend)
        end do
        do bin = 0, nk-1
            signed_bin = bin
            if (.not. real_input .and. bin > (nf-1)/2) signed_bin = bin-nf
            out%frequencies(bin+1) = real(signed_bin,dp)*rate/real(nf,dp)
        end do
        out%power = 0.0_dp
        if (combine == 'adaptive') then
            allocate(pw(nt,nk),q(nt,nk),dof(nk))
            do c = 1, nc
                variance = sum(abs(samples(:,c))**2)/real(n,dp)/rate
                do t = 1, nt
                    transforms = (0.0_dp,0.0_dp)
                    transforms(1:n) = samples(:,c)*windows%tapers(:,t)
                    call fft_state%execute(transforms,fft_status)
                    pw(t,:) = abs(transforms(1:nk))**2/ &
                        (rate*sum(windows%tapers(:,t)**2))
                end do
                call multitaper_weights(pw,windows%ratios,variance,combine,q,dof)
                do t = 1, nt
                    out%power(:,c) = out%power(:,c) + q(t,:)*pw(t,:)
                end do
            end do
            if (real_input) then
                if (mod(nf,2) == 0) then
                    if (nk > 2) out%power(2:nk-1,:) = 2.0_dp*out%power(2:nk-1,:)
                else
                    if (nk > 1) out%power(2:nk,:) = 2.0_dp*out%power(2:nk,:)
                end if
            end if
            out%status = spectral_status_success
            return
        end if
        weight_sum = 0.0_dp
        do t = 1, nt
            wt = 1.0_dp
            if (combine == 'eigen') wt = windows%ratios(t)
            if (wt <= 0.0_dp) cycle
            factor = wt/(rate*sum(windows%tapers(:,t)**2))
            do c = 1, nc
                transforms = (0.0_dp,0.0_dp)
                transforms(1:n) = samples(:,c)*windows%tapers(:,t)
                call fft_state%execute(transforms,fft_status)
                out%power(:,c) = out%power(:,c) + factor*abs(transforms(1:nk))**2
            end do
            weight_sum = weight_sum + wt
        end do
        if (weight_sum <= 0.0_dp) return
        out%power = out%power/weight_sum
        if (real_input) then
            if (mod(nf,2) == 0) then
                if (nk > 2) out%power(2:nk-1,:) = 2.0_dp*out%power(2:nk-1,:)
            else
                if (nk > 1) out%power(2:nk,:) = 2.0_dp*out%power(2:nk,:)
            end if
        end if
        out%status = spectral_status_success
    end function multitaper_channels_core

    function multitaper_psd_real_channels(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode) result(out)
        real(dp), intent(in) :: x(:,:) !! matrix of real samples
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        integer, intent(in), optional :: axis !! time dimension, default 1
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_channels_result) :: out
        integer :: ax
        ax = 1
        if (present(axis)) ax = axis
        select case (ax)
        case (1)
            out = multitaper_channels_core(cmplx(x,0.0_dp,dp),.true.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
        case (2)
            out = multitaper_channels_core(cmplx(transpose(x),0.0_dp,dp), &
                .true.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
        case default
            allocate(out%frequencies(0),out%power(0,0))
        end select
    end function multitaper_psd_real_channels

    function multitaper_psd_complex_channels(x,nw,kmax,fs,nfft,detrend,weighting,axis,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:,:) !! matrix of complex samples
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        integer, intent(in), optional :: axis !! time dimension, default 1
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_channels_result) :: out
        integer :: ax
        ax = 1
        if (present(axis)) ax = axis
        select case (ax)
        case (1)
            out = multitaper_channels_core(x,.false.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
        case (2)
            out = multitaper_channels_core(transpose(x),.false.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
        case default
            allocate(out%frequencies(0),out%power(0,0))
        end select
    end function multitaper_psd_complex_channels

    function multitaper_psd_complex(x,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! complex samples
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove mean
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(spectral_result) :: out
        type(multitaper_matrix_result) :: matrix
        matrix = multitaper_matrix_complex(reshape(x,[size(x),1]),nw,kmax,fs,nfft,detrend,weighting, &
            detrend_mode=detrend_mode)
        out%frequencies = matrix%frequencies
        if (matrix%status == spectral_status_success) then
            out%power = real(matrix%power(:,1,1),dp)
        else
            allocate(out%power(0))
        end if
        out%status = matrix%status
    end function multitaper_psd_complex

    ! A frequency-local median estimates the broad colored-noise background
    ! without requiring a parametric AR/ARMA model. Reflect endpoints for
    ! one-sided spectra and wrap for unshifted full FFT spectra.
    subroutine median_noise_background(power,one_sided,median_bins,background)
        real(dp), intent(in) :: power(:) !! nonnegative PSD samples
        logical, intent(in) :: one_sided !! reflect or periodic frequency axis
        integer, intent(in) :: median_bins !! odd median window length
        real(dp), intent(out) :: background(:) !! median-filtered PSD
        real(dp), allocatable :: sorted(:)
        real(dp) :: value
        integer :: n, half, i, j, k, loc, reflected

        n = size(power)
        half = median_bins/2
        allocate(sorted(median_bins))
        do i = 1, n
            do j = 1, median_bins
                loc = i+j-half-1
                if (one_sided) then
                    reflected = modulo(loc-1,2*n)
                    if (reflected >= n) reflected = 2*n-1-reflected
                    loc = reflected+1
                else
                    loc = modulo(loc-1,n)+1
                end if
                value = power(loc)
                k = j-1
                do while (k >= 1)
                    if (sorted(k) <= value) exit
                    sorted(k+1) = sorted(k)
                    k = k-1
                end do
                sorted(k+1) = value
            end do
            background(i) = sorted(half+1)
        end do
    end subroutine median_noise_background

    function noise_floor_core(psd,one_sided,median_bins) result(out)
        type(spectral_result), intent(in) :: psd !! multitaper PSD and frequencies
        logical, intent(in) :: one_sided !! real-input frequency convention
        integer, intent(in) :: median_bins !! odd frequency-window size, at least 3
        type(multitaper_noise_floor_result) :: out
        integer :: n

        allocate(out%frequencies(0),out%power(0),out%background(0),out%whitened(0))
        if (psd%status /= spectral_status_success) return
        n = size(psd%power)
        if (median_bins < 3 .or. median_bins > n) return
        if (mod(median_bins,2) /= 1) return
        if (.not. all(ieee_is_finite(psd%power))) return
        if (any(psd%power < 0.0_dp)) return
        out%frequencies = psd%frequencies
        out%power = psd%power
        deallocate(out%background,out%whitened)
        allocate(out%background(n),out%whitened(n))
        call median_noise_background(psd%power,one_sided,median_bins,out%background)
        call normalize_noise_power(psd%power,out%background,out%whitened)
        out%status = spectral_status_success
    end function noise_floor_core

    function multitaper_noise_floor_real(x,nw,median_bins,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! real-valued time-domain samples
        real(dp), intent(in) :: nw !! DPSS time-half-bandwidth product
        integer, intent(in) :: median_bins !! odd local median width, 3 to number of bins
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling frequency, default 1 Hz
        integer, intent(in), optional :: nfft !! FFT length, at least size(x)
        logical, intent(in), optional :: detrend !! subtract mean, default true
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_noise_floor_result) :: out
        type(spectral_result) :: psd

        ! Forwarding an absent optional dummy is permitted in Fortran 2008.
        psd = multitaper_psd_real(x,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
        out = noise_floor_core(psd,.true.,median_bins)
    end function multitaper_noise_floor_real

    function multitaper_noise_floor_complex(x,nw,median_bins,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! complex-valued time-domain samples
        real(dp), intent(in) :: nw !! DPSS time-half-bandwidth product
        integer, intent(in) :: median_bins !! odd local median width, 3 to number of bins
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling frequency, default 1 Hz
        integer, intent(in), optional :: nfft !! FFT length, at least size(x)
        logical, intent(in), optional :: detrend !! subtract mean, default true
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_noise_floor_result) :: out
        type(spectral_result) :: psd

        psd = multitaper_psd_complex(x,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
        out = noise_floor_core(psd,.false.,median_bins)
    end function multitaper_noise_floor_complex

    ! Apply a local-median frequency background to each complete window and
    ! channel. For real input, mirror at the ends of the one-sided spectrum;
    ! for complex input, wrap the two-sided FFT ordering periodically.
    function noise_floor_spectrogram_core(raw,one_sided,median_bins) result(out)
        type(multitaper_spectrogram_channels_result), intent(in) :: raw !! input PSD cube
        logical, intent(in) :: one_sided !! one-sided or two-sided FFT
        integer, intent(in) :: median_bins !! odd frequency-window width >= 3
        type(multitaper_noise_floor_spectrogram_channels_result) :: out
        integer :: nf, nt, nc, frame, channel

        allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
            out%background(0,0,0),out%whitened(0,0,0))
        if (raw%status /= spectral_status_success) return
        nf = size(raw%power,1)
        nt = size(raw%power,2)
        nc = size(raw%power,3)
        if (median_bins < 3 .or. median_bins > nf .or. mod(median_bins,2) /= 1) return
        if (.not. all(ieee_is_finite(raw%power))) return
        if (any(raw%power < 0.0_dp)) return
        if (size(raw%frequencies) /= nf .or. size(raw%times) /= nt) return
        out%frequencies = raw%frequencies
        out%times = raw%times
        out%power = raw%power
        deallocate(out%background,out%whitened)
        allocate(out%background(nf,nt,nc),out%whitened(nf,nt,nc))
        do channel = 1, nc
            do frame = 1, nt
                call median_noise_background(raw%power(:,frame,channel),one_sided, &
                    median_bins,out%background(:,frame,channel))
                call normalize_noise_power(raw%power(:,frame,channel), &
                    out%background(:,frame,channel),out%whitened(:,frame,channel))
            end do
        end do
        out%status = spectral_status_success
    end function noise_floor_spectrogram_core

    ! Avoid Inf/NaN in the zero-background convention, matching the
    ! whole-record estimator. This is a descriptive ratio, not a p-value.
    subroutine normalize_noise_power(power,background,whitened)
        real(dp), intent(in) :: power(:) !! raw power >= 0
        real(dp), intent(in) :: background(:) !! median background >= 0
        real(dp), intent(out) :: whitened(:) !! normalized power
        integer :: i

        do i = 1, size(power)
            if (power(i) == 0.0_dp) then
                whitened(i) = 0.0_dp
            else if (background(i) <= 0.0_dp) then
                whitened(i) = huge(1.0_dp)
            else if (background(i) < 1.0_dp .and. &
                power(i) > huge(1.0_dp)*background(i)) then
                whitened(i) = huge(1.0_dp)
            else
                whitened(i) = power(i)/background(i)
            end if
        end do
    end subroutine normalize_noise_power

    function noise_floor_spectrogram_single(multi) result(out)
        type(multitaper_noise_floor_spectrogram_channels_result), intent(in) :: multi !! one channel
        type(multitaper_noise_floor_spectrogram_result) :: out

        out%frequencies = multi%frequencies
        out%times = multi%times
        if (multi%status == spectral_status_success) then
            out%power = multi%power(:,:,1)
            out%background = multi%background(:,:,1)
            out%whitened = multi%whitened(:,:,1)
        else
            allocate(out%power(0,0),out%background(0,0),out%whitened(0,0))
        end if
        out%status = multi%status
    end function noise_floor_spectrogram_single

    function multitaper_noise_floor_spectrogram_real(x,nw,median_bins,nperseg,noverlap,kmax, &
        fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! real signal
        real(dp), intent(in) :: nw !! time-half-bandwidth product
        integer, intent(in) :: median_bins !! odd local median width
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap in samples
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean, default true
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_noise_floor_spectrogram_result) :: out
        type(multitaper_spectrogram_result) :: spec
        type(multitaper_spectrogram_channels_result) :: multi

        spec = multitaper_spectrogram_real(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        multi%status = spec%status
        multi%frequencies = spec%frequencies
        multi%times = spec%times
        multi%power = reshape(spec%power,[size(spec%power,1),size(spec%power,2),1])
        out = noise_floor_spectrogram_single(noise_floor_spectrogram_core(multi,.true.,median_bins))
    end function multitaper_noise_floor_spectrogram_real

    function multitaper_noise_floor_spectrogram_complex(x,nw,median_bins,nperseg,noverlap,kmax, &
        fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! complex signal
        real(dp), intent(in) :: nw !! time-half-bandwidth product
        integer, intent(in) :: median_bins !! odd local median width
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap in samples
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean, default true
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_noise_floor_spectrogram_result) :: out
        type(multitaper_spectrogram_result) :: spec
        type(multitaper_spectrogram_channels_result) :: multi

        spec = multitaper_spectrogram_complex(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        multi%status = spec%status
        multi%frequencies = spec%frequencies
        multi%times = spec%times
        multi%power = reshape(spec%power,[size(spec%power,1),size(spec%power,2),1])
        out = noise_floor_spectrogram_single(noise_floor_spectrogram_core(multi,.false.,median_bins))
    end function multitaper_noise_floor_spectrogram_complex

    function multitaper_noise_floor_spectrogram_real_channels(x,nw,median_bins,nperseg,noverlap,kmax, &
        fs,nfft,detrend,weighting,axis,detrend_mode) result(out)
        real(dp), intent(in) :: x(:,:) !! real signal matrix
        real(dp), intent(in) :: nw !! time-half-bandwidth product
        integer, intent(in) :: median_bins !! odd local median width
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap in samples
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! one-based sample dimension, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_noise_floor_spectrogram_channels_result) :: out
        type(multitaper_spectrogram_channels_result) :: spec

        spec = multitaper_spectrogram_real_channels(x,nw,nperseg,noverlap,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode)
        out = noise_floor_spectrogram_core(spec,.true.,median_bins)
    end function multitaper_noise_floor_spectrogram_real_channels

    function multitaper_noise_floor_spectrogram_complex_channels(x,nw,median_bins,nperseg,noverlap,kmax, &
        fs,nfft,detrend,weighting,axis,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:,:) !! complex signal matrix
        real(dp), intent(in) :: nw !! time-half-bandwidth product
        integer, intent(in) :: median_bins !! odd local median width
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap in samples
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract mean
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        integer, intent(in), optional :: axis !! one-based sample dimension, 1 or 2
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_noise_floor_spectrogram_channels_result) :: out
        type(multitaper_spectrogram_channels_result) :: spec

        spec = multitaper_spectrogram_complex_channels(x,nw,nperseg,noverlap,kmax,fs, &
            nfft,detrend,weighting,axis,detrend_mode)
        out = noise_floor_spectrogram_core(spec,.false.,median_bins)
    end function multitaper_noise_floor_spectrogram_complex_channels

    function multitaper_csd_core(x,y,real_inputs,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        logical, intent(in) :: real_inputs !! both signals real
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! number of tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(cross_spectral_result) :: out
        type(multitaper_matrix_result) :: matrix
        complex(dp), allocatable :: data(:,:)

        allocate(out%frequencies(0),out%power(0))
        if (size(x) /= size(y)) return
        allocate(data(size(x),2))
        data(:,1) = x
        data(:,2) = y
        matrix = multitaper_matrix_core(data,real_inputs,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
        if (matrix%status == spectral_status_success) then
            out%frequencies = matrix%frequencies
            out%power = matrix%power(:,1,2)
        end if
        out%status = matrix%status
    end function multitaper_csd_core

    function multitaper_csd_real(x,y,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! first real signal
        real(dp), intent(in) :: y(:) !! second real signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(cross_spectral_result) :: out
        out = multitaper_csd_core(cmplx(x,0.0_dp,dp),cmplx(y,0.0_dp,dp), &
            .true.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_csd_real

    function multitaper_csd_complex(x,y,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first complex signal
        complex(dp), intent(in) :: y(:) !! second complex signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(cross_spectral_result) :: out
        out = multitaper_csd_core(x,y,.false.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_csd_complex

    function multitaper_csd_real_complex(x,y,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! first real signal
        complex(dp), intent(in) :: y(:) !! second complex signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(cross_spectral_result) :: out
        out = multitaper_csd_core(cmplx(x,0.0_dp,dp),y,.false.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_csd_real_complex

    function multitaper_csd_complex_real(x,y,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first complex signal
        real(dp), intent(in) :: y(:) !! second real signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! tapers
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(cross_spectral_result) :: out
        out = multitaper_csd_core(x,cmplx(y,0.0_dp,dp),.false.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_csd_complex_real

    function multitaper_coherence_core(x,y,real_inputs,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        logical, intent(in) :: real_inputs !! both signals real
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(spectral_result) :: out
        type(multitaper_matrix_result) :: matrix
        complex(dp), allocatable :: data(:,:)
        real(dp) :: denom
        integer :: k

        allocate(out%frequencies(0),out%power(0))
        if (size(x) /= size(y)) return
        allocate(data(size(x),2))
        data(:,1) = x
        data(:,2) = y
        matrix = multitaper_matrix_core(data,real_inputs,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
        if (matrix%status /= spectral_status_success) return
        out%frequencies = matrix%frequencies
        deallocate(out%power)
        allocate(out%power(size(out%frequencies)))
        do k = 1, size(out%power)
            denom = real(matrix%power(k,1,1)*matrix%power(k,2,2),dp)
            out%power(k) = 0.0_dp
            if (denom > 0.0_dp) out%power(k) = min(1.0_dp, &
                abs(matrix%power(k,1,2))**2/denom)
        end do
        out%status = spectral_status_success
    end function multitaper_coherence_core

    function multitaper_coherence_real(x,y,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! first real signal
        real(dp), intent(in) :: y(:) !! second real signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(spectral_result) :: out
        out = multitaper_coherence_core(cmplx(x,0.0_dp,dp),cmplx(y,0.0_dp,dp), &
            .true.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_coherence_real

    function multitaper_coherence_complex(x,y,nw,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first complex signal
        complex(dp), intent(in) :: y(:) !! second complex signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! transform length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal or eigen
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(spectral_result) :: out
        out = multitaper_coherence_core(x,y,.false.,nw,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_coherence_complex


    ! Cached-taper overlap spectrogram: each segment is normalized as a PSD,
    ! not an STFT magnitude, to preserve multitaper_psd units and frequency
    ! bins. Odd and even real FFT lengths have different Nyquist doubling.
    function multitaper_spectrogram_core(x,real_input,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:,:) !! samples, channels
        logical, intent(in) :: real_input !! one-sided spectrum for real data
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! number of samples in each window
        integer, intent(in), optional :: noverlap !! overlapping samples; default nperseg/2
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! sampling frequency, default 1
        integer, intent(in), optional :: nfft !! transform length, at least nperseg
        logical, intent(in), optional :: detrend !! subtract segment mean, default true
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_spectrogram_channels_result) :: out

        type(dpss_result) :: windows
        type(fft_plan) :: fft_state
        integer :: fft_status
        complex(dp), allocatable :: segment(:), transform(:)
        real(dp), allocatable :: periodograms(:,:), weights(:,:), dof(:)
        real(dp) :: rate, variance, norm_taper
        integer :: n, nc, nf, nk, nt, overlap, hop, nframe, frame, ch, taper, bin, j, start, signed_bin
        character(:), allocatable :: mode, trend

        allocate(out%frequencies(0),out%times(0),out%power(0,0,0))
        n = size(x,1)
        nc = size(x,2)
        if (nperseg < 2 .or. nperseg > n .or. nc < 1) return
        if (.not. ieee_is_finite(nw)) return
        if (nw <= 0.0_dp .or. nw >= 0.5_dp*real(nperseg,dp)) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp) return
        nf = nperseg
        if (present(nfft)) nf = nfft
        if (nf < nperseg) return
        overlap = nperseg/2
        if (present(noverlap)) overlap = noverlap
        if (overlap < 0 .or. overlap >= nperseg) return
        hop = nperseg-overlap
        nt = min(nperseg,max(1,int(2.0_dp*nw)-1))
        if (present(kmax)) nt = kmax
        if (nt < 1 .or. nt > nperseg) return
        mode = 'equal'
        if (present(weighting)) mode = trim(weighting)
        if (mode /= 'equal' .and. mode /= 'eigen' .and. mode /= 'adaptive') return
        trend = 'constant'
        if (present(detrend)) then
            if (.not. detrend) trend = 'none'
        end if
        if (present(detrend_mode)) trend = trim(detrend_mode)
        if (trend /= 'none' .and. trend /= 'constant' .and. trend /= 'linear') return

        ! DPSS eigensystem is computed just once for all windows and channels.
        windows = dpss(nperseg,nw,kmax=nt,sym=.false.,norm='l2')
        if (windows%status /= dpss_status_success) return
        call fft_state%initialize(nf,fft_status)
        if (fft_status /= spectral_status_success) return
        nframe = 1+(n-nperseg)/hop
        nk = nf
        if (real_input) nk = nf/2+1
        deallocate(out%frequencies,out%times,out%power)
        allocate(out%frequencies(nk),out%times(nframe),out%power(nk,nframe,nc))
        allocate(segment(nperseg),transform(nf),periodograms(nt,nk),weights(nt,nk),dof(nk))
        do bin = 0, nk-1
            signed_bin = bin
            if (.not. real_input .and. bin > (nf-1)/2) signed_bin = bin-nf
            out%frequencies(bin+1) = real(signed_bin,dp)*rate/real(nf,dp)
        end do
        do frame = 1, nframe
            start = 1+(frame-1)*hop
            ! Match SciPy's legacy spectrogram sample-time convention.
            out%times(frame) = real(start-1+nperseg/2,dp)/rate
            do ch = 1, nc
                segment = x(start:start+nperseg-1,ch)
                call detrend_multitaper_segment(segment,trend)
                variance = sum(abs(segment)**2)/real(nperseg,dp)/rate
                do taper = 1, nt
                    transform = (0.0_dp,0.0_dp)
                    transform(1:nperseg) = segment*windows%tapers(:,taper)
                    call fft_state%execute(transform,fft_status)
                    norm_taper = rate*sum(windows%tapers(:,taper)**2)
                    periodograms(taper,:) = abs(transform(1:nk))**2/norm_taper
                end do
                call multitaper_weights(periodograms,windows%ratios,variance,mode,weights,dof)
                do bin = 1, nk
                    out%power(bin,frame,ch) = sum(weights(:,bin)*periodograms(:,bin))
                end do
                if (real_input) then
                    j = nk
                    if (mod(nf,2) == 0) j = nk-1
                    if (j >= 2) out%power(2:j,frame,ch) = 2.0_dp*out%power(2:j,frame,ch)
                end if
            end do
        end do
        out%status = spectral_status_success
    end function multitaper_spectrogram_core

    function multitaper_spectrogram_real(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! real signal
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! samples per second
        integer, intent(in), optional :: nfft !! Fourier transform length
        logical, intent(in), optional :: detrend !! subtract mean per window
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_spectrogram_result) :: out
        type(multitaper_spectrogram_channels_result) :: multi

        multi = multitaper_spectrogram_core(reshape(cmplx(x,0.0_dp,dp),[size(x),1]), &
            .true.,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        out%frequencies = multi%frequencies
        out%times = multi%times
        if (multi%status == spectral_status_success) then
            out%power = multi%power(:,:,1)
        else
            allocate(out%power(0,0))
        end if
        out%status = multi%status
    end function multitaper_spectrogram_real

    function multitaper_spectrogram_complex(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! complex signal
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! samples per second
        integer, intent(in), optional :: nfft !! Fourier transform length
        logical, intent(in), optional :: detrend !! subtract mean per window
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_spectrogram_result) :: out
        type(multitaper_spectrogram_channels_result) :: multi

        multi = multitaper_spectrogram_core(reshape(x,[size(x),1]), &
            .false.,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        out%frequencies = multi%frequencies
        out%times = multi%times
        if (multi%status == spectral_status_success) then
            out%power = multi%power(:,:,1)
        else
            allocate(out%power(0,0))
        end if
        out%status = multi%status
    end function multitaper_spectrogram_complex

    function multitaper_spectrogram_real_channels(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting, &
        axis,detrend_mode) result(out)
        real(dp), intent(in) :: x(:,:) !! matrix of channels
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! samples per second
        integer, intent(in), optional :: nfft !! Fourier transform length
        logical, intent(in), optional :: detrend !! subtract mean per window
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        integer, intent(in), optional :: axis !! one-based sample dimension, default 1
        type(multitaper_spectrogram_channels_result) :: out
        integer :: sample_axis

        sample_axis = 1
        if (present(axis)) sample_axis = axis
        if (sample_axis == 1) then
            out = multitaper_spectrogram_core(cmplx(x,0.0_dp,dp), &
                .true.,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        else if (sample_axis == 2) then
            out = multitaper_spectrogram_core(cmplx(transpose(x),0.0_dp,dp), &
                .true.,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        else
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0))
        end if
    end function multitaper_spectrogram_real_channels

    function multitaper_spectrogram_complex_channels(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend, &
        weighting,axis,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:,:) !! matrix of channels
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in) :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! number of DPSS tapers
        real(dp), intent(in), optional :: fs !! samples per second
        integer, intent(in), optional :: nfft !! Fourier transform length
        logical, intent(in), optional :: detrend !! subtract mean per window
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        integer, intent(in), optional :: axis !! one-based sample dimension, default 1
        type(multitaper_spectrogram_channels_result) :: out
        integer :: sample_axis

        sample_axis = 1
        if (present(axis)) sample_axis = axis
        if (sample_axis == 1) then
            out = multitaper_spectrogram_core(x,.false.,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        else if (sample_axis == 2) then
            out = multitaper_spectrogram_core(transpose(x), &
                .false.,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        else
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0))
        end if
    end function multitaper_spectrogram_complex_channels

    ! Thomson's harmonic F-test, fitted independently at each DFT frequency.
    ! Under the usual locally Gaussian smooth-spectrum null, the statistic
    ! follows F(2,2*K-2) approximately at interior frequencies. P-values are
    ! pointwise, not corrected for multiple frequency comparisons. The fitted
    ! complex amplitude equals half the physical sinusoid amplitude for a
    ! real signal at an isolated positive interior Fourier frequency.
    function multitaper_harmonic_core(x,real_input,nw,kmax,fs,nfft,detrend,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! input signal
        logical, intent(in) :: real_input !! use one-sided frequencies
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in), optional :: kmax !! number of DPSS tapers (at least 2)
        real(dp), intent(in), optional :: fs !! sampling rate
        integer, intent(in), optional :: nfft !! FFT length, at least size(x)
        logical, intent(in), optional :: detrend !! mean removal, default false
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_harmonic_result) :: out

        type(dpss_result) :: windows
        type(fft_plan) :: fft_state
        integer :: fft_status
        complex(dp), allocatable :: samples(:), transform(:), tapered(:,:)
        complex(dp) :: numerator, fitted
        real(dp), allocatable :: response(:)
        real(dp) :: rate, response_norm, rss, explained, fval, pval
        integer :: n, nf, nk, nt, taper, bin, signed_bin
        character(:), allocatable :: trend

        allocate(out%frequencies(0),out%f_statistic(0),out%pvalue(0),out%amplitude(0))
        n = size(x)
        if (n < 2 .or. .not. ieee_is_finite(nw)) return
        if (nw <= 0.0_dp .or. nw >= 0.5_dp*real(n,dp)) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp) return
        nf = n
        if (present(nfft)) nf = nfft
        if (nf < n) return
        nt = min(n,max(2,int(2.0_dp*nw)-1))
        if (present(kmax)) nt = kmax
        if (nt < 2 .or. nt > n) return
        trend = 'none'
        if (present(detrend)) then
            if (detrend) then
                trend = 'constant'
            else
                trend = 'none'
            end if
        end if
        if (present(detrend_mode)) trend = trim(detrend_mode)
        if (trend /= 'none' .and. trend /= 'constant' .and. trend /= 'linear') return
        windows = dpss(n,nw,kmax=nt,sym=.false.,norm='l2')
        if (windows%status /= dpss_status_success) return
        call fft_state%initialize(nf,fft_status)
        if (fft_status /= spectral_status_success) return

        nk = nf
        if (real_input) nk = nf/2+1
        deallocate(out%frequencies,out%f_statistic,out%pvalue,out%amplitude)
        allocate(out%frequencies(nk),out%f_statistic(nk),out%pvalue(nk),out%amplitude(nk))
        allocate(samples(n),transform(nf),tapered(nt,nk),response(nt))
        samples = x
        call detrend_multitaper_segment(samples,trend)
        response = sum(windows%tapers,dim=1)
        response_norm = sum(response**2)
        if (response_norm <= tiny(1.0_dp)) return
        do taper = 1, nt
            transform = (0.0_dp,0.0_dp)
            transform(1:n) = samples*windows%tapers(:,taper)
            call fft_state%execute(transform,fft_status)
            tapered(taper,:) = transform(1:nk)
        end do
        do bin = 0, nk-1
            signed_bin = bin
            if (.not. real_input .and. bin > (nf-1)/2) signed_bin = bin-nf
            out%frequencies(bin+1) = real(signed_bin,dp)*rate/real(nf,dp)
            numerator = sum(response*tapered(:,bin+1))
            fitted = numerator/response_norm
            rss = sum(abs(tapered(:,bin+1)-fitted*response)**2)
            explained = abs(numerator)**2/response_norm
            fval = 0.0_dp
            pval = 1.0_dp
            if (explained > 0.0_dp) then
                if (rss <= tiny(1.0_dp)) then
                    fval = huge(1.0_dp)
                    pval = 0.0_dp
                else
                    fval = real(nt-1,dp)*explained/rss
                    pval = exp(-real(nt-1,dp)*log(1.0_dp+fval/real(nt-1,dp)))
                end if
            end if
            out%f_statistic(bin+1) = fval
            out%pvalue(bin+1) = pval
            out%amplitude(bin+1) = fitted
        end do
        out%status = spectral_status_success
    end function multitaper_harmonic_core

    function multitaper_harmonic_ftest_real(x,nw,kmax,fs,nfft,detrend,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! real input signal
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in), optional :: kmax !! number of tapers, at least 2
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract mean, default false
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_harmonic_result) :: out
        out = multitaper_harmonic_core(cmplx(x,0.0_dp,dp),.true.,nw,kmax,fs,nfft,detrend,detrend_mode)
    end function multitaper_harmonic_ftest_real

    function multitaper_harmonic_ftest_complex(x,nw,kmax,fs,nfft,detrend,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! complex input signal
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in), optional :: kmax !! number of tapers, at least 2
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract mean, default false
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_harmonic_result) :: out
        out = multitaper_harmonic_core(x,.false.,nw,kmax,fs,nfft,detrend,detrend_mode)
    end function multitaper_harmonic_ftest_complex
    ! Stable bottom-up merge sort of frequency indices by increasing p-value.
    ! O(m log m), unlike insertion sorting long zero-padded FFT grids.
    subroutine sort_pvalue_indices(p,indices)
        real(dp), intent(in) :: p(:) !! pointwise probabilities
        integer, intent(inout) :: indices(:) !! selected one-based bins
        integer, allocatable :: work(:)
        integer :: n, width, left, mid, right, i, j, k

        n = size(indices)
        allocate(work(n))
        width = 1
        do while (width < n)
            left = 1
            do while (left <= n)
                mid = left + min(width,n-left+1) - 1
                right = mid + min(width,n-mid)
                i = left
                j = mid+1
                do k = left, right
                    if (i > mid) then
                        work(k) = indices(j)
                        j = j+1
                    else if (j > right) then
                        work(k) = indices(i)
                        i = i+1
                    else if (p(indices(i)) <= p(indices(j))) then
                        work(k) = indices(i)
                        i = i+1
                    else
                        work(k) = indices(j)
                        j = j+1
                    end if
                end do
                left = right+1
            end do
            indices = work
            if (width > n/2) exit
            width = 2*width
        end do
    end subroutine sort_pvalue_indices

    ! Correct the existing Thomson pointwise F-test, without altering the
    ! regression statistic or assuming independent Fourier bins.
    function harmonic_scan_core(raw,real_input,nfft,method,alpha,include_edges,test_bins) result(out)
        type(multitaper_harmonic_result), intent(in) :: raw !! unadjusted line tests
        logical, intent(in) :: real_input !! one-sided Fourier frequencies
        integer, intent(in) :: nfft !! underlying FFT length
        character(*), intent(in), optional :: method !! none, bonferroni, holm, bh, by
        real(dp), intent(in), optional :: alpha !! significance level, default 0.05
        logical, intent(in), optional :: include_edges !! test DC/Nyquist, default false
        integer, intent(in), optional :: test_bins(:) !! explicitly selected one-based bins
        type(multitaper_harmonic_scan_result) :: out

        integer, allocatable :: indices(:)
        logical, allocatable :: selected(:)
        character(:), allocatable :: correction
        real(dp) :: cutoff, multiplier, running, candidate
        integer :: n, m, i, j, bin
        logical :: include_boundary

        allocate(out%frequencies(0),out%f_statistic(0),out%pvalue(0), &
            out%adjusted_pvalue(0),out%amplitude(0),out%significant(0))
        if (raw%status /= spectral_status_success) return
        cutoff = 0.05_dp
        if (present(alpha)) cutoff = alpha
        if (.not. ieee_is_finite(cutoff)) return
        if (cutoff <= 0.0_dp .or. cutoff >= 1.0_dp) return
        correction = 'bh'
        if (present(method)) correction = trim(method)
        select case (correction)
        case ('none','bonferroni','holm','bh','by')
        case default
            return
        end select
        n = size(raw%pvalue)
        if (n < 1 .or. nfft < 1) return
        if (any(.not. ieee_is_finite(raw%pvalue))) return
        if (any(raw%pvalue < 0.0_dp .or. raw%pvalue > 1.0_dp)) return
        include_boundary = .false.
        if (present(include_edges)) include_boundary = include_edges
        allocate(selected(n))
        selected = .true.
        if (present(test_bins)) then
            selected = .false.
            do i = 1, size(test_bins)
                bin = test_bins(i)
                if (bin < 1 .or. bin > n) return
                if (selected(bin)) return
                selected(bin) = .true.
            end do
        else if (.not. include_boundary) then
            selected(1) = .false.  ! DC violates the interior-bin F null model
            if (real_input .and. mod(nfft,2) == 0) selected(nfft/2+1) = .false.
        end if
        m = count(selected)
        allocate(indices(m))
        j = 0
        do i = 1, n
            if (.not. selected(i)) cycle
            j = j+1
            indices(j) = i
        end do
        deallocate(out%frequencies,out%f_statistic,out%pvalue, &
            out%adjusted_pvalue,out%amplitude,out%significant)
        out%frequencies = raw%frequencies
        out%f_statistic = raw%f_statistic
        out%pvalue = raw%pvalue
        out%amplitude = raw%amplitude
        allocate(out%adjusted_pvalue(n),out%significant(n))
        out%adjusted_pvalue = 1.0_dp
        out%significant = .false.
        out%num_tests = m
        if (m == 0) then
            out%status = spectral_status_success
            return
        end if
        select case (correction)
        case ('none')
            do j = 1, m
                i = indices(j)
                out%adjusted_pvalue(i) = raw%pvalue(i)
            end do
        case ('bonferroni')
            do j = 1, m
                i = indices(j)
                out%adjusted_pvalue(i) = min(1.0_dp,real(m,dp)*raw%pvalue(i))
            end do
        case ('holm')
            ! Holm step-down FWER correction; valid under arbitrary dependence.
            call sort_pvalue_indices(raw%pvalue,indices)
            running = 0.0_dp
            do j = 1, m
                i = indices(j)
                candidate = min(1.0_dp,real(m-j+1,dp)*raw%pvalue(i))
                running = max(running,candidate)
                out%adjusted_pvalue(i) = running
            end do
        case ('bh','by')
            call sort_pvalue_indices(raw%pvalue,indices)
            multiplier = 1.0_dp
            if (correction == 'by') then
                multiplier = 0.0_dp
                do j = 1, m
                    multiplier = multiplier+1.0_dp/real(j,dp)
                end do
            end if
            running = 1.0_dp
            do j = m, 1, -1
                i = indices(j)
                candidate = min(1.0_dp,raw%pvalue(i)*(real(m,dp)/real(j,dp))*multiplier)
                running = min(running,candidate)
                out%adjusted_pvalue(i) = running
            end do
        end select
        out%significant = selected .and. out%adjusted_pvalue <= cutoff
        out%status = spectral_status_success
    end function harmonic_scan_core

    function multitaper_harmonic_scan_real(x,nw,kmax,fs,nfft,detrend,method,alpha,include_edges,test_bins,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! sampled real signal
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in), optional :: kmax !! DPSS taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT size
        logical, intent(in), optional :: detrend !! subtract mean
        character(*), intent(in), optional :: method !! none, bonferroni, holm, bh, by
        real(dp), intent(in), optional :: alpha !! target error rate
        logical, intent(in), optional :: include_edges !! include DC/Nyquist tests
        integer, intent(in), optional :: test_bins(:) !! explicit one-based bins
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_harmonic_scan_result) :: out
        type(multitaper_harmonic_result) :: raw
        integer :: nf

        nf = size(x)
        if (present(nfft)) nf = nfft
        raw = multitaper_harmonic_ftest_real(x,nw,kmax,fs,nfft,detrend,detrend_mode)
        out = harmonic_scan_core(raw,.true.,nf,method,alpha,include_edges,test_bins)
    end function multitaper_harmonic_scan_real

    function multitaper_harmonic_scan_complex(x,nw,kmax,fs,nfft,detrend,method,alpha,include_edges, &
        test_bins,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! sampled complex signal
        real(dp), intent(in) :: nw !! time-bandwidth product
        integer, intent(in), optional :: kmax !! DPSS taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT size
        logical, intent(in), optional :: detrend !! subtract mean
        character(*), intent(in), optional :: method !! none, bonferroni, holm, bh, by
        real(dp), intent(in), optional :: alpha !! target error rate
        logical, intent(in), optional :: include_edges !! include DC test
        integer, intent(in), optional :: test_bins(:) !! explicit one-based bins
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear; overrides detrend
        type(multitaper_harmonic_scan_result) :: out
        type(multitaper_harmonic_result) :: raw
        integer :: nf

        nf = size(x)
        if (present(nfft)) nf = nfft
        raw = multitaper_harmonic_ftest_complex(x,nw,kmax,fs,nfft,detrend,detrend_mode)
        out = harmonic_scan_core(raw,.false.,nf,method,alpha,include_edges,test_bins)
    end function multitaper_harmonic_scan_complex

    ! Multitaper cross-spectral matrices, resolved in time. Both channels
    ! independently choose adaptive weights. Their geometric mean gives a
    ! Cauchy-Schwarz-consistent cross spectrum and coherence in [0,1].
    ! Frequencies follow the unshifted DFT order, with the usual one-sided
    ! doubling for real inputs except DC and the even-length Nyquist bin.
    function multitaper_matrix_spectrogram_core(x,real_input,nw,nperseg,noverlap,kmax,fs,nfft,detrend, &
        weighting,detrend_mode,cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! samples by channels
        logical, intent(in) :: real_input !! both signals are real
        real(dp), intent(in) :: nw !! time-half-bandwidth product
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlapping samples
        integer, intent(in), optional :: kmax !! DPSS taper count
        real(dp), intent(in), optional :: fs !! sampling rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! subtract segment means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        logical, intent(in), optional :: cv_shrinkage !! evaluate LOO spectra for shrinkage
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        type(dpss_result) :: tapers
        type(fft_plan) :: fft_state
        integer :: fft_status
        complex(dp), allocatable :: segment(:), transform(:), spectra(:,:,:)
        real(dp), allocatable :: periodograms(:,:), weights(:,:,:), dof(:), variances(:)
        real(dp) :: rate, variance, normalization
        logical :: cross_validate
        integer :: n, nc, nf, nk, nt, overlap, hop, nframe, frame, ch, c, d, taper, bin, start, signed_bin, j
        character(:), allocatable :: mode, trend

        allocate(out%frequencies(0),out%times(0),out%power(0,0,0,0))
        n = size(x,1)
        nc = size(x,2)
        if (nperseg < 2 .or. nperseg > n .or. nc < 1) return
        if (.not. ieee_is_finite(nw)) return
        if (nw <= 0.0_dp .or. nw >= 0.5_dp*real(nperseg,dp)) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp) return
        nf = nperseg
        if (present(nfft)) nf = nfft
        if (nf < nperseg) return
        overlap = nperseg/2
        if (present(noverlap)) overlap = noverlap
        if (overlap < 0 .or. overlap >= nperseg) return
        hop = nperseg-overlap
        nt = min(nperseg,max(1,int(2.0_dp*nw)-1))
        if (present(kmax)) nt = kmax
        if (nt < 1 .or. nt > nperseg) return
        mode = 'equal'
        if (present(weighting)) mode = trim(weighting)
        if (mode /= 'equal' .and. mode /= 'eigen' .and. mode /= 'adaptive') return
        cross_validate = .false.
        if (present(cv_shrinkage)) cross_validate = cv_shrinkage
        if (cross_validate .and. (nt < 3 .or. nc < 2)) return
        trend = 'constant'
        if (present(detrend)) then
            if (.not. detrend) trend = 'none'
        end if
        if (present(detrend_mode)) trend = trim(detrend_mode)
        if (trend /= 'none' .and. trend /= 'constant' .and. trend /= 'linear') return

        ! Calculate the tapers once, not for every frame or every channel.
        tapers = dpss(nperseg,nw,kmax=nt,sym=.false.,norm='l2')
        if (tapers%status /= dpss_status_success) return
        call fft_state%initialize(nf,fft_status)
        if (fft_status /= spectral_status_success) return
        nframe = 1+(n-nperseg)/hop
        nk = nf
        if (real_input) nk = nf/2+1
        deallocate(out%frequencies,out%times,out%power)
        allocate(out%frequencies(nk),out%times(nframe),out%power(nk,nframe,nc,nc))
        allocate(segment(nperseg),transform(nf),spectra(nk,nt,nc))
        allocate(periodograms(nt,nk),weights(nt,nk,nc),dof(nk),variances(nc))
        if (cross_validate) allocate(out%cv_shrinkage(nk,nframe))
        do bin = 0, nk-1
            signed_bin = bin
            if (.not. real_input .and. bin > (nf-1)/2) signed_bin = bin-nf
            out%frequencies(bin+1) = real(signed_bin,dp)*rate/real(nf,dp)
        end do
        out%power = (0.0_dp,0.0_dp)
        do frame = 1, nframe
            start = 1+(frame-1)*hop
            out%times(frame) = real(start-1+nperseg/2,dp)/rate
            do ch = 1, nc
                segment = x(start:start+nperseg-1,ch)
                call detrend_multitaper_segment(segment,trend)
                variance = sum(abs(segment)**2)/real(nperseg,dp)/rate
                variances(ch) = variance
                do taper = 1, nt
                    transform = (0.0_dp,0.0_dp)
                    transform(1:nperseg) = segment*tapers%tapers(:,taper)
                    call fft_state%execute(transform,fft_status)
                    normalization = sqrt(rate*sum(tapers%tapers(:,taper)**2))
                    spectra(:,taper,ch) = transform(1:nk)/normalization
                    periodograms(taper,:) = abs(spectra(:,taper,ch))**2
                end do
                call multitaper_weights(periodograms,tapers%ratios,variance,mode,weights(:,:,ch),dof)
            end do
            if (cross_validate) then
                call multitaper_cv_diagonal(spectra,tapers%ratios,variances,mode,out%cv_shrinkage(:,frame))
            end if
            do d = 1, nc
                do c = 1, nc
                    do taper = 1, nt
                        out%power(:,frame,c,d) = out%power(:,frame,c,d) + &
                            sqrt(weights(taper,:,c)*weights(taper,:,d))* &
                            conjg(spectra(:,taper,c))*spectra(:,taper,d)
                    end do
                end do
            end do
            if (real_input) then
                j = nk
                if (mod(nf,2) == 0) j = nk-1
                if (j >= 2) out%power(2:j,frame,:,:) = 2.0_dp*out%power(2:j,frame,:,:)
            end if
        end do
        if (mode == 'equal') out%effective_tapers = real(nt,dp)
        if (mode == 'eigen') then
            out%effective_tapers = multitaper_effective_tapers(max(tapers%ratios,0.0_dp))
            ! Match the time-window equal-weight fallback for zero eigen ratios.
            if (out%effective_tapers == 0.0_dp) out%effective_tapers = real(nt,dp)
        end if
        out%status = spectral_status_success
    end function multitaper_matrix_spectrogram_core

    function multitaper_matrix_spectrogram_real(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting, &
        axis,detrend_mode,cv_shrinkage) result(out)
        real(dp), intent(in) :: x(:,:) !! input real channels
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! segment size
        integer, intent(in), optional :: noverlap !! overlap count
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! mean removal
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        integer, intent(in), optional :: axis !! sample axis (1 or 2)
        logical, intent(in), optional :: cv_shrinkage !! optional LOO shrinkage estimate
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        integer :: ax
        ax = 1
        if (present(axis)) ax = axis
        if (ax == 1) then
            out = multitaper_matrix_spectrogram_core(cmplx(x,0.0_dp,dp),.true., &
                nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode,cv_shrinkage)
        else if (ax == 2) then
            out = multitaper_matrix_spectrogram_core(cmplx(transpose(x),0.0_dp,dp),.true., &
                nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode,cv_shrinkage)
        else
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0,0))
        end if
    end function multitaper_matrix_spectrogram_real

    function multitaper_matrix_spectrogram_complex(x,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting, &
        axis,detrend_mode,cv_shrinkage) result(out)
        complex(dp), intent(in) :: x(:,:) !! input complex channels
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! segment size
        integer, intent(in), optional :: noverlap !! overlap count
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! mean removal
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        integer, intent(in), optional :: axis !! sample axis (1 or 2)
        logical, intent(in), optional :: cv_shrinkage !! optional LOO shrinkage estimate
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        integer :: ax
        ax = 1
        if (present(axis)) ax = axis
        if (ax == 1) then
            out = multitaper_matrix_spectrogram_core(x,.false., &
                nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode,cv_shrinkage)
        else if (ax == 2) then
            out = multitaper_matrix_spectrogram_core(transpose(x),.false., &
                nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode,cv_shrinkage)
        else
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0,0))
        end if
    end function multitaper_matrix_spectrogram_complex

    function multitaper_cross_spectrogram_core(x,y,real_inputs,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) &
        result(out)
        complex(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        logical, intent(in) :: real_inputs !! both input signals real
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_cross_spectrogram_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: matrix
        complex(dp), allocatable :: samples(:,:)

        allocate(out%frequencies(0),out%times(0),out%power(0,0),out%auto_x(0,0),out%auto_y(0,0))
        if (size(x) /= size(y)) return
        allocate(samples(size(x),2))
        samples(:,1) = x
        samples(:,2) = y
        matrix = multitaper_matrix_spectrogram_core(samples,real_inputs,nw,nperseg,noverlap,kmax, &
            fs,nfft,detrend,weighting,detrend_mode)
        if (matrix%status /= spectral_status_success) return
        out%frequencies = matrix%frequencies
        out%times = matrix%times
        out%power = matrix%power(:,:,1,2)
        out%auto_x = real(matrix%power(:,:,1,1),dp)
        out%auto_y = real(matrix%power(:,:,2,2),dp)
        out%status = spectral_status_success
    end function multitaper_cross_spectrogram_core

    function multitaper_tf_coherence_core(x,y,real_inputs,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) &
        result(out)
        complex(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        logical, intent(in) :: real_inputs !! both input signals real
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_tf_coherence_result) :: out
        type(multitaper_cross_spectrogram_result) :: cross
        real(dp), allocatable :: denom(:,:)

        allocate(out%frequencies(0),out%times(0),out%coherence(0,0))
        cross = multitaper_cross_spectrogram_core(x,y,real_inputs,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
        if (cross%status /= spectral_status_success) return
        out%frequencies = cross%frequencies
        out%times = cross%times
        allocate(denom(size(cross%auto_x,1),size(cross%auto_x,2)))
        denom = max(cross%auto_x,0.0_dp)*max(cross%auto_y,0.0_dp)
        deallocate(out%coherence)
        allocate(out%coherence(size(denom,1),size(denom,2)))
        out%coherence = 0.0_dp
        where (denom > tiny(1.0_dp))
            out%coherence = min(1.0_dp,abs(cross%power)**2/max(denom,tiny(1.0_dp)))
        end where
        out%status = spectral_status_success
    end function multitaper_tf_coherence_core

    function multitaper_cross_spectrogram_real(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! first signal
        real(dp), intent(in) :: y(:) !! second signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_cross_spectrogram_result) :: out
        out = multitaper_cross_spectrogram_core(cmplx(x,0.0_dp,dp),cmplx(y,0.0_dp,dp),.true.,nw, &
            nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_cross_spectrogram_real

    function multitaper_cross_spectrogram_complex(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_cross_spectrogram_result) :: out
        out = multitaper_cross_spectrogram_core(x,y,.false.,nw, &
            nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_cross_spectrogram_complex

    function multitaper_cross_spectrogram_real_complex(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend, &
        weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_cross_spectrogram_result) :: out
        out = multitaper_cross_spectrogram_core(cmplx(x,0.0_dp,dp),y,.false.,nw, &
            nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_cross_spectrogram_real_complex

    function multitaper_cross_spectrogram_complex_real(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend, &
        weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first signal
        real(dp), intent(in) :: y(:) !! second signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_cross_spectrogram_result) :: out
        out = multitaper_cross_spectrogram_core(x,cmplx(y,0.0_dp,dp),.false.,nw, &
            nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_cross_spectrogram_complex_real

    function multitaper_tf_coherence_real(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! first signal
        real(dp), intent(in) :: y(:) !! second signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_tf_coherence_result) :: out
        out = multitaper_tf_coherence_core(cmplx(x,0.0_dp,dp),cmplx(y,0.0_dp,dp),.true.,nw, &
            nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_tf_coherence_real

    function multitaper_tf_coherence_complex(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_tf_coherence_result) :: out
        out = multitaper_tf_coherence_core(x,y,.false.,nw, &
            nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_tf_coherence_complex

    function multitaper_tf_coherence_real_complex(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        real(dp), intent(in) :: x(:) !! first signal
        complex(dp), intent(in) :: y(:) !! second signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_tf_coherence_result) :: out
        out = multitaper_tf_coherence_core(cmplx(x,0.0_dp,dp),y,.false.,nw, &
            nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_tf_coherence_real_complex

    function multitaper_tf_coherence_complex_real(x,y,nw,nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode) result(out)
        complex(dp), intent(in) :: x(:) !! first signal
        real(dp), intent(in) :: y(:) !! second signal
        real(dp), intent(in) :: nw !! time-half-bandwidth
        integer, intent(in) :: nperseg !! window length
        integer, intent(in), optional :: noverlap !! overlap samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sample rate
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove means
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear
        type(multitaper_tf_coherence_result) :: out
        out = multitaper_tf_coherence_core(x,cmplx(y,0.0_dp,dp),.false.,nw, &
            nperseg,noverlap,kmax,fs,nfft,detrend,weighting,detrend_mode)
    end function multitaper_tf_coherence_complex_real


    ! Initialize a streaming spectral estimator once; DPSS tapers and work
    ! arrays are retained across calls. Call reset() to restart the signal
    ! without recomputing the tapers. The optional complex_input flag chooses
    ! two-sided output; otherwise real one-sided spectra are produced.
    subroutine multitaper_stream_initialize(self,nw,nperseg,nchannels,status,noverlap,kmax,fs,nfft, &
        detrend,weighting,complex_input, &
        detrend_mode)
        class(multitaper_stream), intent(inout) :: self
        real(dp), intent(in) :: nw !! time-half-bandwidth product
        integer, intent(in) :: nperseg !! window length
        integer, intent(in) :: nchannels !! number of simultaneous input channels
        integer, intent(out) :: status !! spectral_status_success or invalid_input
        integer, intent(in), optional :: noverlap !! overlapping samples
        integer, intent(in), optional :: kmax !! taper count
        real(dp), intent(in), optional :: fs !! sampling frequency
        integer, intent(in), optional :: nfft !! FFT length
        logical, intent(in), optional :: detrend !! remove mean per channel and window
        character(*), intent(in), optional :: weighting !! equal/eigen/adaptive
        logical, intent(in), optional :: complex_input !! use two-sided spectrum
        character(*), intent(in), optional :: detrend_mode !! none/constant/linear per complete window
        type(dpss_result) :: w
        integer :: overlap, taper, bin, signed_bin
        character(:), allocatable :: mode
        real(dp) :: rate

        self%ready = .false.
        call self%reset()
        if (allocated(self%weighting)) deallocate(self%weighting)
        if (allocated(self%trend)) deallocate(self%trend)
        if (allocated(self%frequencies)) deallocate(self%frequencies)
        if (allocated(self%tapers)) deallocate(self%tapers)
        if (allocated(self%ratios)) deallocate(self%ratios)
        if (allocated(self%normalizers)) deallocate(self%normalizers)
        if (allocated(self%periodograms)) deallocate(self%periodograms)
        if (allocated(self%weights)) deallocate(self%weights)
        if (allocated(self%dof)) deallocate(self%dof)
        if (allocated(self%window)) deallocate(self%window)
        if (allocated(self%segment)) deallocate(self%segment)
        if (allocated(self%transform)) deallocate(self%transform)
        if (allocated(self%spectra)) deallocate(self%spectra)
        if (allocated(self%summed_spectra)) deallocate(self%summed_spectra)
        call self%fft_state%clear()
        status = spectral_status_invalid_input
        if (nperseg < 2 .or. nchannels < 1) return
        if (.not. ieee_is_finite(nw)) return
        if (nw <= 0.0_dp .or. nw >= 0.5_dp*real(nperseg,dp)) return
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp) return
        overlap = nperseg/2
        if (present(noverlap)) overlap = noverlap
        if (overlap < 0 .or. overlap >= nperseg) return
        self%nt = min(nperseg,max(1,int(2.0_dp*nw)-1))
        if (present(kmax)) self%nt = kmax
        if (self%nt < 1 .or. self%nt > nperseg) return
        self%nfft = nperseg
        if (present(nfft)) self%nfft = nfft
        if (self%nfft < nperseg) return
        mode = 'equal'
        if (present(weighting)) mode = trim(weighting)
        if (mode /= 'equal' .and. mode /= 'eigen' .and. mode /= 'adaptive') return
        self%real_input = .true.
        if (present(complex_input)) self%real_input = .not. complex_input
        self%trend = 'constant'
        if (present(detrend)) then
            if (.not. detrend) self%trend = 'none'
        end if
        if (present(detrend_mode)) self%trend = trim(detrend_mode)
        if (self%trend /= 'none' .and. self%trend /= 'constant' .and. self%trend /= 'linear') return
        w = dpss(nperseg,nw,kmax=self%nt,sym=.false.,norm='l2')
        if (w%status /= dpss_status_success) return
        call self%fft_state%initialize(self%nfft,status)
        if (status /= spectral_status_success) return

        self%nperseg = nperseg
        self%nc = nchannels
        self%rate = rate
        self%hop = nperseg-overlap
        self%weighting = mode
        self%nk = self%nfft
        if (self%real_input) self%nk = self%nfft/2+1
        self%tapers = w%tapers
        self%ratios = w%ratios
        allocate(self%normalizers(self%nt),self%frequencies(self%nk))
        allocate(self%window(nperseg,nchannels),self%segment(nperseg),self%transform(self%nfft))
        allocate(self%spectra(self%nk,self%nt,nchannels),self%periodograms(self%nt,self%nk))
        allocate(self%weights(self%nt,self%nk,nchannels),self%dof(self%nk))
        allocate(self%summed_spectra(self%nk,nchannels,nchannels))
        do taper = 1, self%nt
            self%normalizers(taper) = sqrt(rate*sum(self%tapers(:,taper)**2))
        end do
        do bin = 0, self%nk-1
            signed_bin = bin
            if (.not. self%real_input .and. bin > (self%nfft-1)/2) signed_bin = bin-self%nfft
            self%frequencies(bin+1) = real(signed_bin,dp)*rate/real(self%nfft,dp)
        end do
        self%ready = .true.
        call self%reset()
        status = spectral_status_success
    end subroutine multitaper_stream_initialize

    subroutine multitaper_stream_reset(self)
        class(multitaper_stream), intent(inout) :: self
        self%buffered = 0
        self%samples_seen = 0_int64
        self%frames_seen = 0_int64
        if (allocated(self%window)) self%window = (0.0_dp,0.0_dp)
        if (allocated(self%summed_spectra)) self%summed_spectra = (0.0_dp,0.0_dp)
    end subroutine multitaper_stream_reset

    integer(int64) function multitaper_stream_sample_count(self) result(n)
        class(multitaper_stream), intent(in) :: self
        n = self%samples_seen
    end function multitaper_stream_sample_count

    integer(int64) function multitaper_stream_window_count(self) result(n)
        class(multitaper_stream), intent(in) :: self
        n = self%frames_seen
    end function multitaper_stream_window_count

    ! Allocate exactly the new complete windows. No result is retained by the
    ! state after the caller receives it, including the time-frequency cube.
    subroutine multitaper_stream_output_shell(self,nchunk,out)
        class(multitaper_stream), intent(in) :: self
        integer, intent(in) :: nchunk !! new samples on each channel
        type(multitaper_spectral_matrix_spectrogram_result), intent(out) :: out
        integer :: nframes, f
        allocate(out%frequencies(0),out%times(0),out%power(0,0,0,0))
        if (.not. self%ready) return
        nframes = 0
        if (nchunk >= self%nperseg-self%buffered) then
            nframes = 1+(nchunk-self%nperseg+self%buffered)/self%hop
        end if
        deallocate(out%frequencies,out%times,out%power)
        allocate(out%frequencies(self%nk),out%times(nframes))
        allocate(out%power(self%nk,nframes,self%nc,self%nc))
        out%frequencies = self%frequencies
        do f = 1, nframes
            out%times(f) = (real(self%frames_seen+int(f, int64)-1_int64,dp)*real(self%hop,dp) &
                + real(self%nperseg/2,dp))/self%rate
        end do
        out%power = (0.0_dp,0.0_dp)
        out%status = spectral_status_success
    end subroutine multitaper_stream_output_shell

    ! This is the same window estimator used by the batch cross-spectrogram:
    ! independently weighted channel spectra and geometrically combined
    ! cross-channel taper weights, including one-sided real FFT doubling.
    subroutine multitaper_stream_compute_window(self,matrix)
        class(multitaper_stream), intent(inout) :: self
        complex(dp), intent(out) :: matrix(:,:,:) !! frequency, channel, channel
        real(dp) :: variance
        integer :: ch, taper, c, d, j, fft_status

        do ch = 1, self%nc
            self%segment = self%window(:,ch)
            call detrend_multitaper_segment(self%segment,self%trend)
            variance = sum(abs(self%segment)**2)/real(self%nperseg,dp)/self%rate
            do taper = 1, self%nt
                self%transform = (0.0_dp,0.0_dp)
                self%transform(1:self%nperseg) = self%segment*self%tapers(:,taper)
                call self%fft_state%execute(self%transform,fft_status)
                if (fft_status /= spectral_status_success) then
                    matrix = (0.0_dp,0.0_dp)
                    return
                end if
                self%spectra(:,taper,ch) = self%transform(1:self%nk)/self%normalizers(taper)
                self%periodograms(taper,:) = abs(self%spectra(:,taper,ch))**2
            end do
            call multitaper_weights(self%periodograms,self%ratios,variance,self%weighting, &
                self%weights(:,:,ch),self%dof)
        end do
        matrix = (0.0_dp,0.0_dp)
        do d = 1, self%nc
            do c = 1, self%nc
                do taper = 1, self%nt
                    matrix(:,c,d) = matrix(:,c,d) + &
                        sqrt(self%weights(taper,:,c)*self%weights(taper,:,d))* &
                        conjg(self%spectra(:,taper,c))*self%spectra(:,taper,d)
                end do
            end do
        end do
        if (self%real_input) then
            j = self%nk
            if (mod(self%nfft,2) == 0) j = self%nk-1
            if (j >= 2) matrix(2:j,:,:) = 2.0_dp*matrix(2:j,:,:)
        end if
    end subroutine multitaper_stream_compute_window

    subroutine multitaper_stream_emit(self,out,iframe)
        class(multitaper_stream), intent(inout) :: self
        type(multitaper_spectral_matrix_spectrogram_result), intent(inout) :: out
        integer, intent(inout) :: iframe
        integer :: overlap
        if (self%buffered /= self%nperseg) return
        iframe = iframe+1
        call self%compute_window(out%power(:,iframe,:,:))
        self%summed_spectra = self%summed_spectra+out%power(:,iframe,:,:)
        self%frames_seen = self%frames_seen+1_int64
        overlap = self%nperseg-self%hop
        if (overlap > 0) self%window(1:overlap,:) = self%window(self%hop+1:self%nperseg,:)
        self%buffered = overlap
    end subroutine multitaper_stream_emit

    function multitaper_stream_push_real(self,x) result(out)
        class(multitaper_stream), intent(inout) :: self
        real(dp), intent(in) :: x(:) !! next real single-channel data chunk
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        integer :: i, frame
        call self%output_shell(size(x),out)
        if (out%status /= spectral_status_success) return
        if (self%nc /= 1 .or. .not. self%real_input .or. .not. all(ieee_is_finite(x))) then
            out%status = spectral_status_invalid_input
            return
        end if
        frame = 0
        do i = 1, size(x)
            self%buffered = self%buffered+1
            self%window(self%buffered,1) = cmplx(x(i),0.0_dp,dp)
            self%samples_seen = self%samples_seen+1_int64
            call self%emit(out,frame)
        end do
    end function multitaper_stream_push_real

    function multitaper_stream_push_complex(self,x) result(out)
        class(multitaper_stream), intent(inout) :: self
        complex(dp), intent(in) :: x(:) !! next complex single-channel data chunk
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        integer :: i, frame
        call self%output_shell(size(x),out)
        if (out%status /= spectral_status_success) return
        if (self%nc /= 1 .or. self%real_input .or. &
            .not. all(ieee_is_finite(real(x,dp))) .or. .not. all(ieee_is_finite(aimag(x)))) then
            out%status = spectral_status_invalid_input
            return
        end if
        frame = 0
        do i = 1, size(x)
            self%buffered = self%buffered+1
            self%window(self%buffered,1) = x(i)
            self%samples_seen = self%samples_seen+1_int64
            call self%emit(out,frame)
        end do
    end function multitaper_stream_push_complex

    function multitaper_stream_push_real_channels(self,x,axis) result(out)
        class(multitaper_stream), intent(inout) :: self
        real(dp), intent(in) :: x(:,:) !! new samples, channels
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        integer, intent(in), optional :: axis !! one-based sample dimension: 1 (default) or 2
        integer :: i, frame, ax, count, channels
        ax = 1
        if (present(axis)) ax = axis
        count = 0
        channels = 0
        if (ax == 1) then
            count = size(x,1)
            channels = size(x,2)
        else if (ax == 2) then
            count = size(x,2)
            channels = size(x,1)
        end if
        call self%output_shell(count,out)
        if (ax /= 1 .and. ax /= 2) then
            out%status = spectral_status_invalid_input
            return
        end if
        if (out%status /= spectral_status_success) return
        if (channels /= self%nc .or. .not. self%real_input .or. .not. all(ieee_is_finite(x))) then
            out%status = spectral_status_invalid_input
            return
        end if
        frame = 0
        do i = 1, count
            self%buffered = self%buffered+1
            if (ax == 1) then
                self%window(self%buffered,:) = cmplx(x(i,:),0.0_dp,dp)
            else
                self%window(self%buffered,:) = cmplx(x(:,i),0.0_dp,dp)
            end if
            self%samples_seen = self%samples_seen+1_int64
            call self%emit(out,frame)
        end do
    end function multitaper_stream_push_real_channels

    function multitaper_stream_push_complex_channels(self,x,axis) result(out)
        class(multitaper_stream), intent(inout) :: self
        complex(dp), intent(in) :: x(:,:) !! new complex samples, channels
        type(multitaper_spectral_matrix_spectrogram_result) :: out
        integer, intent(in), optional :: axis !! one-based sample dimension: 1 (default) or 2
        integer :: i, frame, ax, count, channels
        ax = 1
        if (present(axis)) ax = axis
        count = 0
        channels = 0
        if (ax == 1) then
            count = size(x,1)
            channels = size(x,2)
        else if (ax == 2) then
            count = size(x,2)
            channels = size(x,1)
        end if
        call self%output_shell(count,out)
        if (ax /= 1 .and. ax /= 2) then
            out%status = spectral_status_invalid_input
            return
        end if
        if (out%status /= spectral_status_success) return
        if (channels /= self%nc .or. self%real_input .or. &
            .not. all(ieee_is_finite(real(x,dp))) .or. .not. all(ieee_is_finite(aimag(x)))) then
            out%status = spectral_status_invalid_input
            return
        end if
        frame = 0
        do i = 1, count
            self%buffered = self%buffered+1
            if (ax == 1) then
                self%window(self%buffered,:) = x(i,:)
            else
                self%window(self%buffered,:) = x(:,i)
            end if
            self%samples_seen = self%samples_seen+1_int64
            call self%emit(out,frame)
        end do
    end function multitaper_stream_push_complex_channels

    ! Running mean cross-spectral matrix over *complete* windows only. It is
    ! kept in the state; the caller need not keep any previously emitted bins.
    function multitaper_stream_mean_matrix(self) result(out)
        class(multitaper_stream), intent(in) :: self
        type(multitaper_matrix_result) :: out
        allocate(out%frequencies(0),out%power(0,0,0))
        if (.not. self%ready .or. self%frames_seen == 0_int64) return
        deallocate(out%frequencies,out%power)
        allocate(out%frequencies(self%nk),out%power(self%nk,self%nc,self%nc))
        out%frequencies = self%frequencies
        out%power = self%summed_spectra/real(self%frames_seen,dp)
        out%status = spectral_status_success
    end function multitaper_stream_mean_matrix

    ! Convert only the diagonal of the current streaming cross-spectral cube
    ! to channel PSDs; no history of prior frames is required.
    function multitaper_stream_floor_from_matrix(self,mat,median_bins) result(out)
        class(multitaper_stream), intent(in) :: self
        type(multitaper_spectral_matrix_spectrogram_result), intent(in) :: mat !! emitted frames
        integer, intent(in) :: median_bins !! odd median-frequency width
        type(multitaper_noise_floor_spectrogram_channels_result) :: out
        type(multitaper_spectrogram_channels_result) :: psd
        integer :: ch

        psd%status = mat%status
        psd%frequencies = mat%frequencies
        psd%times = mat%times
        allocate(psd%power(size(mat%power,1),size(mat%power,2),size(mat%power,3)))
        if (mat%status == spectral_status_success) then
            do ch = 1, size(psd%power,3)
                psd%power(:,:,ch) = max(0.0_dp,real(mat%power(:,:,ch,ch),dp))
            end do
        end if
        out = noise_floor_spectrogram_core(psd,self%real_input,median_bins)
    end function multitaper_stream_floor_from_matrix

    function multitaper_stream_floor_real(self,x,median_bins) result(out)
        class(multitaper_stream), intent(inout) :: self
        real(dp), intent(in) :: x(:) !! next real samples
        integer, intent(in) :: median_bins !! odd local median width
        type(multitaper_noise_floor_spectrogram_channels_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: mat

        if (.not. self%ready) then
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
                out%background(0,0,0),out%whitened(0,0,0))
            return
        end if
        if (median_bins < 3 .or. median_bins > self%nk .or. mod(median_bins,2) /= 1) then
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
                out%background(0,0,0),out%whitened(0,0,0))
            return
        end if
        mat = self%push_real(x)
        out = multitaper_stream_floor_from_matrix(self,mat,median_bins)
    end function multitaper_stream_floor_real

    function multitaper_stream_floor_complex(self,x,median_bins) result(out)
        class(multitaper_stream), intent(inout) :: self
        complex(dp), intent(in) :: x(:) !! next complex samples
        integer, intent(in) :: median_bins !! odd local median width
        type(multitaper_noise_floor_spectrogram_channels_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: mat

        if (.not. self%ready) then
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
                out%background(0,0,0),out%whitened(0,0,0))
            return
        end if
        if (median_bins < 3 .or. median_bins > self%nk .or. mod(median_bins,2) /= 1) then
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
                out%background(0,0,0),out%whitened(0,0,0))
            return
        end if
        mat = self%push_complex(x)
        out = multitaper_stream_floor_from_matrix(self,mat,median_bins)
    end function multitaper_stream_floor_complex

    function multitaper_stream_floor_real_channels(self,x,median_bins,axis) result(out)
        class(multitaper_stream), intent(inout) :: self
        real(dp), intent(in) :: x(:,:) !! new real samples, channels
        integer, intent(in) :: median_bins !! odd local median width
        integer, intent(in), optional :: axis !! sample dimension, 1 or 2
        type(multitaper_noise_floor_spectrogram_channels_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: mat

        if (.not. self%ready) then
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
                out%background(0,0,0),out%whitened(0,0,0))
            return
        end if
        if (median_bins < 3 .or. median_bins > self%nk .or. mod(median_bins,2) /= 1) then
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
                out%background(0,0,0),out%whitened(0,0,0))
            return
        end if
        mat = self%push_real_channels(x,axis)
        out = multitaper_stream_floor_from_matrix(self,mat,median_bins)
    end function multitaper_stream_floor_real_channels

    function multitaper_stream_floor_complex_channels(self,x,median_bins,axis) result(out)
        class(multitaper_stream), intent(inout) :: self
        complex(dp), intent(in) :: x(:,:) !! new complex samples, channels
        integer, intent(in) :: median_bins !! odd local median width
        integer, intent(in), optional :: axis !! sample dimension, 1 or 2
        type(multitaper_noise_floor_spectrogram_channels_result) :: out
        type(multitaper_spectral_matrix_spectrogram_result) :: mat

        if (.not. self%ready) then
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
                out%background(0,0,0),out%whitened(0,0,0))
            return
        end if
        if (median_bins < 3 .or. median_bins > self%nk .or. mod(median_bins,2) /= 1) then
            allocate(out%frequencies(0),out%times(0),out%power(0,0,0), &
                out%background(0,0,0),out%whitened(0,0,0))
            return
        end if
        mat = self%push_complex_channels(x,axis)
        out = multitaper_stream_floor_from_matrix(self,mat,median_bins)
    end function multitaper_stream_floor_complex_channels

    ! Normalize the channel PSDs of the running mean spectral matrix.
    ! The median is applied after the frame mean, not to each frame before it.
    function multitaper_stream_mean_noise_floor(self,median_bins) result(out)
        class(multitaper_stream), intent(in) :: self
        integer, intent(in) :: median_bins !! odd local median frequency width
        type(multitaper_noise_floor_channels_result) :: out
        type(multitaper_matrix_result) :: avg
        integer :: nf, nc, ch

        allocate(out%frequencies(0),out%power(0,0),out%background(0,0),out%whitened(0,0))
        avg = self%mean_matrix()
        if (avg%status /= spectral_status_success) return
        nf = size(avg%power,1)
        nc = size(avg%power,2)
        if (median_bins < 3 .or. median_bins > nf .or. mod(median_bins,2) /= 1) return
        out%frequencies = avg%frequencies
        deallocate(out%power,out%background,out%whitened)
        allocate(out%power(nf,nc),out%background(nf,nc),out%whitened(nf,nc))
        do ch = 1, nc
            out%power(:,ch) = max(0.0_dp,real(avg%power(:,ch,ch),dp))
            call median_noise_background(out%power(:,ch),self%real_input,median_bins,out%background(:,ch))
            call normalize_noise_power(out%power(:,ch),out%background(:,ch),out%whitened(:,ch))
        end do
        out%status = spectral_status_success
    end function multitaper_stream_mean_noise_floor

end module scifort_signal_multitaper
