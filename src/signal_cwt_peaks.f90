! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Ricker continuous-wavelet peak detection. Wavelet/ridge/noise conventions
! are adapted from SciPy 1.17 scipy.signal._peak_finding (BSD-3-Clause).
! SciPy copyright and license are reproduced in THIRD_PARTY_LICENSES.md.
! See CODE_PROVENANCE.md for the Bioinformatics ridge-detection citation.
module scifort_signal_cwt_peaks
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private
    integer, parameter, public :: cwt_peak_success = 0
    integer, parameter, public :: cwt_peak_invalid_input = 1
    public :: ricker_wavelet, cwt_ricker, cwt_wavelet, find_peaks_cwt
    public :: cwt_wavelet_callback
    public :: morlet2_wavelet, cwt_morlet2, cwt_complex_wavelet
    public :: cwt_complex_wavelet_callback
    public :: find_peaks_cwt_complex, find_peaks_cwt_magnitude

    interface cwt_complex_wavelet
        module procedure cwt_complex_wavelet_real
        module procedure cwt_complex_wavelet_complex
    end interface cwt_complex_wavelet

    interface cwt_morlet2
        module procedure cwt_morlet2_real
        module procedure cwt_morlet2_complex
    end interface cwt_morlet2

    interface find_peaks_cwt_complex
        module procedure find_peaks_cwt_complex_real
        module procedure find_peaks_cwt_complex_complex
    end interface find_peaks_cwt_complex

    abstract interface
        function cwt_wavelet_callback(points, width) result(w)
            import dp
            real(dp), intent(in) :: points !! requested number of samples, possibly noninteger
            real(dp), intent(in) :: width !! positive wavelet scale
            real(dp), allocatable :: w(:) !! finite real wavelet, ceiling(points) samples
        end function cwt_wavelet_callback
        function cwt_complex_wavelet_callback(points, width) result(w)
            import dp
            real(dp), intent(in) :: points !! unrounded requested wavelet length
            real(dp), intent(in) :: width !! strictly positive wavelet scale
            complex(dp), allocatable :: w(:) !! complex wavelet, ceiling(points) samples
        end function cwt_complex_wavelet_callback
    end interface
contains
    pure function ricker_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! wavelet length before rounding, as in SciPy
        real(dp), intent(in) :: width !! positive scale parameter
        real(dp), allocatable :: w(:) !! sampled Mexican-hat wavelet
        real(dp) :: x, scale, t
        integer :: k, n

        allocate(w(0))
        if (.not. ieee_is_finite(points) .or. .not. ieee_is_finite(width)) return
        if (points <= 0.0_dp .or. width <= 0.0_dp) return
        if (points > real(huge(n) - 1, dp)) return
        n = ceiling(points)
        deallocate(w)
        allocate(w(n))
        scale = 2.0_dp / (sqrt(3.0_dp * width) * acos(-1.0_dp)**0.25_dp)
        do k = 1, n
            x = (real(k - 1, dp) - (points - 1.0_dp) / 2.0_dp) / width
            t = x*x
            w(k) = scale * (1.0_dp - t) * exp(-0.5_dp * t)
        end do
    end function ricker_wavelet

    function cwt_ricker(x, widths, status) result(cwt)
        real(dp), intent(in) :: x(:) !! finite input signal
        real(dp), intent(in) :: widths(:) !! positive wavelet widths
        integer, intent(out), optional :: status !! zero on success
        real(dp), allocatable :: cwt(:,:) !! wavelet row, signal column

        cwt = cwt_wavelet(x, widths, ricker_wavelet, status)
    end function cwt_ricker

    function cwt_wavelet(x, widths, wavelet, status) result(cwt)
        real(dp), intent(in) :: x(:) !! finite input signal
        real(dp), intent(in) :: widths(:) !! finite positive wavelet widths
        procedure(cwt_wavelet_callback) :: wavelet !! callback accepting (points, width)
        integer, intent(out), optional :: status !! zero on success, one for bad data or kernel
        real(dp), allocatable :: cwt(:,:) !! wavelet row, signal column; empty on error
        real(dp), allocatable :: wave(:)
        real(dp) :: points
        integer :: i, j, k, ns, nw, m, offset, index

        if (present(status)) status = cwt_peak_invalid_input
        allocate(cwt(0,0))
        if (size(x) == 0 .or. size(widths) == 0) return
        if (any(.not. ieee_is_finite(x))) return
        if (any(.not. ieee_is_finite(widths))) return
        if (any(widths <= 0.0_dp)) return
        if (any(widths > real(huge(ns),dp)/10.0_dp)) return
        ns = size(x)
        nw = size(widths)
        deallocate(cwt)
        allocate(cwt(nw,ns))
        cwt = 0.0_dp
        do i = 1, nw
            points = min(10.0_dp * widths(i), real(ns,dp))
            wave = wavelet(points, widths(i))
            ! SciPy calls wavelet(points,width) with possibly fractional points.
            ! The returned kernel must contain the requested number of samples.
            if (size(wave) /= ceiling(points)) then
                deallocate(cwt)
                allocate(cwt(0,0))
                return
            end if
            if (any(.not. ieee_is_finite(wave))) then
                deallocate(cwt)
                allocate(cwt(0,0))
                return
            end if
            m = size(wave)
            offset = (m - 1)/2
            do j = 1, ns
                do k = 1, m
                    index = j + offset + 1 - k
                    if (index < 1 .or. index > ns) cycle
                    ! SciPy convolves with the wavelet reversed.
                    cwt(i,j) = cwt(i,j) + x(index)*wave(m+1-k)
                end do
            end do
        end do
        if (present(status)) status = cwt_peak_success
    end function cwt_wavelet

    ! The complex CWT uses the conjugate time-reversed wavelet, matching
    ! scipy.signal._wavelets._cwt.  This is deliberately a separate API from
    ! the real-valued ridge peak finder: complex coefficients have no natural
    ! ordering for finding row maxima.
    pure function morlet2_wavelet(points, width) result(w)
        real(dp), intent(in) :: points !! wavelet length before rounding (can be fractional)
        real(dp), intent(in) :: width !! positive scale parameter
        complex(dp), allocatable :: w(:) !! complex Morlet wavelet at omega0=5

        w = morlet2_samples(points, width, 5.0_dp)
    end function morlet2_wavelet

    pure function morlet2_samples(points, width, omega0) result(w)
        real(dp), intent(in) :: points !! requested unrounded length
        real(dp), intent(in) :: width !! time-frequency scale (>0)
        real(dp), intent(in) :: omega0 !! dimensionless angular center frequency
        complex(dp), allocatable :: w(:)
        real(dp) :: x, amplitude
        integer :: n, k

        allocate(w(0))
        if (.not. ieee_is_finite(points) .or. .not. ieee_is_finite(width) .or. &
                .not. ieee_is_finite(omega0)) return
        if (points <= 0.0_dp .or. width <= 0.0_dp) return
        if (points > real(huge(n)-1,dp)) return
        n = ceiling(points)
        deallocate(w)
        allocate(w(n))
        amplitude = acos(-1.0_dp)**(-0.25_dp)/sqrt(width)
        do k = 1, n
            x = (real(k-1,dp)-0.5_dp*(points-1.0_dp))/width
            w(k) = amplitude*exp(-0.5_dp*x*x)*cmplx(cos(omega0*x),sin(omega0*x),dp)
        end do
    end function morlet2_samples

    function cwt_complex_wavelet_real(x, widths, wavelet, status) result(cwt)
        real(dp), intent(in) :: x(:) !! finite real signal
        real(dp), intent(in) :: widths(:) !! positive wavelet scales
        procedure(cwt_complex_wavelet_callback) :: wavelet !! complex wavelet callback
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: cwt(:,:) !! scale by position complex coefficients

        call complex_wavelet_kernel(cmplx(x,0.0_dp,dp), widths, cwt, &
            wavelet=wavelet, status=status)
    end function cwt_complex_wavelet_real

    function cwt_complex_wavelet_complex(x, widths, wavelet, status) result(cwt)
        complex(dp), intent(in) :: x(:) !! finite complex signal
        real(dp), intent(in) :: widths(:) !! positive wavelet scales
        procedure(cwt_complex_wavelet_callback) :: wavelet !! complex wavelet callback
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: cwt(:,:) !! scale by position complex coefficients

        call complex_wavelet_kernel(x, widths, cwt, wavelet=wavelet, status=status)
    end function cwt_complex_wavelet_complex

    function cwt_morlet2_real(x, widths, omega0, status) result(cwt)
        real(dp), intent(in) :: x(:) !! finite real signal
        real(dp), intent(in) :: widths(:) !! strictly positive wavelet scales
        real(dp), intent(in), optional :: omega0 !! angular Morlet center frequency, default 5
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: cwt(:,:) !! scale by position complex coefficients

        call complex_wavelet_kernel(cmplx(x,0.0_dp,dp), widths, cwt, &
            omega0=omega0, status=status)
    end function cwt_morlet2_real

    function cwt_morlet2_complex(x, widths, omega0, status) result(cwt)
        complex(dp), intent(in) :: x(:) !! finite complex signal
        real(dp), intent(in) :: widths(:) !! strictly positive wavelet scales
        real(dp), intent(in), optional :: omega0 !! angular Morlet center frequency, default 5
        integer, intent(out), optional :: status !! success or invalid input
        complex(dp), allocatable :: cwt(:,:) !! scale by position complex coefficients

        call complex_wavelet_kernel(x, widths, cwt, omega0=omega0, status=status)
    end function cwt_morlet2_complex

    subroutine complex_wavelet_kernel(x, widths, cwt, wavelet, omega0, status)
        complex(dp), intent(in) :: x(:) !! signal in the complex domain
        real(dp), intent(in) :: widths(:) !! strictly positive wavelet scales
        complex(dp), allocatable, intent(out) :: cwt(:,:) !! complex CWT (empty on error)
        procedure(cwt_complex_wavelet_callback), optional :: wavelet !! sampled custom wavelet
        real(dp), intent(in), optional :: omega0 !! Morlet angular center frequency
        integer, intent(out), optional :: status !! success or invalid input

        complex(dp), allocatable :: wave(:), computed(:,:)
        real(dp) :: points, frequency
        integer :: i, j, k, ns, nw, m, offset, index

        allocate(cwt(0,0))
        if (present(status)) status = cwt_peak_invalid_input
        ns = size(x)
        nw = size(widths)
        if (ns < 1 .or. nw < 1) return
        if (.not. all(ieee_is_finite(real(x,dp))) .or. &
                .not. all(ieee_is_finite(aimag(x)))) return
        if (.not. all(ieee_is_finite(widths))) return
        if (any(widths <= 0.0_dp)) return
        if (any(widths > real(huge(ns),dp)/10.0_dp)) return
        if (present(wavelet) .and. present(omega0)) return
        frequency = 5.0_dp
        if (present(omega0)) frequency = omega0
        if (.not. ieee_is_finite(frequency)) return

        allocate(computed(nw,ns))
        computed = (0.0_dp,0.0_dp)
        do i = 1, nw
            points = min(10.0_dp*widths(i),real(ns,dp))
            if (present(wavelet)) then
                wave = wavelet(points,widths(i))
            else
                wave = morlet2_samples(points,widths(i),frequency)
            end if
            if (size(wave) /= ceiling(points)) return
            if (.not. all(ieee_is_finite(real(wave,dp))) .or. &
                    .not. all(ieee_is_finite(aimag(wave)))) return
            m = size(wave)
            offset = (m-1)/2
            do j = 1, ns
                do k = 1, m
                    index = j+offset+1-k
                    if (index < 1 .or. index > ns) cycle
                    computed(i,j) = computed(i,j) + x(index)*conjg(wave(m+1-k))
                end do
            end do
        end do
        if (.not. all(ieee_is_finite(real(computed,dp))) .or. &
                .not. all(ieee_is_finite(aimag(computed)))) return
        call move_alloc(computed,cwt)
        if (present(status)) status = cwt_peak_success
    end subroutine complex_wavelet_kernel

    pure function percentile(x, p) result(y)
        real(dp), intent(in) :: x(:) !! finite vector with positive size
        real(dp), intent(in) :: p !! percentile on [0,100]
        real(dp) :: y, q, v
        real(dp) :: sorted(size(x))
        integer :: i, j, lo, hi
        sorted = x
        do i = 2, size(sorted)
            v = sorted(i)
            j = i-1
            do while (j >= 1)
                if (sorted(j) <= v) exit
                sorted(j+1) = sorted(j)
                j = j-1
            end do
            sorted(j+1) = v
        end do
        q = 1.0_dp + (real(size(x)-1,dp)*p/100.0_dp)
        lo = max(1,min(size(x),int(q)))
        hi = min(size(x),lo+1)
        y = sorted(lo) + (q - real(lo,dp))*(sorted(hi)-sorted(lo))
    end function percentile

    function find_peaks_cwt(x, widths, max_distances, gap_thresh, min_length, &
            min_snr, noise_perc, window_size, status, wavelet) result(peaks)
        real(dp), intent(in) :: x(:) !! finite input signal
        real(dp), intent(in) :: widths(:) !! positive wavelet scales
        real(dp), intent(in), optional :: max_distances(:) !! maximum ridge shifts
        real(dp), intent(in), optional :: gap_thresh !! maximum unconnected rows
        integer, intent(in), optional :: min_length !! minimum points in a ridge
        real(dp), intent(in), optional :: min_snr !! minimum coefficient/noise ratio
        real(dp), intent(in), optional :: noise_perc !! percentile of first-scale row
        integer, intent(in), optional :: window_size !! percentile neighborhood size
        integer, intent(out), optional :: status !! success or invalid input
        procedure(cwt_wavelet_callback), optional :: wavelet !! custom real wavelet
        integer, allocatable :: peaks(:) !! sorted one-based peak positions
        real(dp), allocatable :: cwt(:,:)
        integer :: code

        if (present(status)) status = cwt_peak_invalid_input
        allocate(peaks(0))
        if (size(x) == 0 .or. size(widths) == 0) return
        if (present(wavelet)) then
            cwt = cwt_wavelet(x, widths, wavelet, code)
        else
            cwt = cwt_ricker(x, widths, code)
        end if
        if (code /= cwt_peak_success) return
        peaks = ridge_peaks_real(cwt, widths, max_distances, gap_thresh, &
            min_length, min_snr, noise_perc, window_size, status)
    end function find_peaks_cwt

    ! For complex wavelets, maxima and percentile noise are defined by the
    ! modulus of the complex CWT, NOT the real part or a signed coefficient.
    ! This deliberate extension is distinct from SciPy's Ricker peak finder.
    function find_peaks_cwt_magnitude(coefficients, widths, max_distances, &
            gap_thresh, min_length, min_snr, noise_perc, window_size, status) result(peaks)
        complex(dp), intent(in) :: coefficients(:,:) !! scale x sample complex CWT
        real(dp), intent(in) :: widths(:) !! positive scales, one per matrix row
        real(dp), intent(in), optional :: max_distances(:) !! maximum ridge shifts
        real(dp), intent(in), optional :: gap_thresh !! maximum missing rows
        integer, intent(in), optional :: min_length !! minimum ridge length
        real(dp), intent(in), optional :: min_snr !! minimum modulus/noise ratio
        real(dp), intent(in), optional :: noise_perc !! local modulus percentile
        integer, intent(in), optional :: window_size !! noise neighborhood length
        integer, intent(out), optional :: status !! success or invalid input
        integer, allocatable :: peaks(:) !! sorted one-based peak positions
        real(dp), allocatable :: magnitudes(:,:)

        if (present(status)) status = cwt_peak_invalid_input
        allocate(peaks(0))
        if (size(coefficients,1) /= size(widths)) return
        if (size(coefficients,1) == 0 .or. size(coefficients,2) == 0) return
        if (.not. all(ieee_is_finite(real(coefficients,dp))) .or. &
                .not. all(ieee_is_finite(aimag(coefficients)))) return
        magnitudes = abs(coefficients)
        if (.not. all(ieee_is_finite(magnitudes))) return
        peaks = ridge_peaks_real(magnitudes, widths, max_distances, gap_thresh, &
            min_length, min_snr, noise_perc, window_size, status)
    end function find_peaks_cwt_magnitude

    function find_peaks_cwt_complex_real(x, widths, max_distances, &
            gap_thresh, min_length, min_snr, noise_perc, window_size, status, &
            wavelet, omega0) result(peaks)
        real(dp), intent(in) :: x(:) !! real-valued input signal
        real(dp), intent(in) :: widths(:) !! positive complex-wavelet scales
        real(dp), intent(in), optional :: max_distances(:) !! maximum ridge shifts
        real(dp), intent(in), optional :: gap_thresh !! maximum missing rows
        integer, intent(in), optional :: min_length !! minimum ridge length
        real(dp), intent(in), optional :: min_snr !! minimum modulus/noise ratio
        real(dp), intent(in), optional :: noise_perc !! local modulus percentile
        integer, intent(in), optional :: window_size !! noise neighborhood length
        integer, intent(out), optional :: status !! success or invalid input
        procedure(cwt_complex_wavelet_callback), optional :: wavelet !! custom complex kernel
        real(dp), intent(in), optional :: omega0 !! default Morlet angular frequency
        integer, allocatable :: peaks(:) !! sorted one-based peak locations

        call complex_peak_kernel(cmplx(x,0.0_dp,dp), widths, peaks, &
            max_distances, gap_thresh, min_length, min_snr, noise_perc, &
            window_size, status, wavelet, omega0)
    end function find_peaks_cwt_complex_real

    function find_peaks_cwt_complex_complex(x, widths, max_distances, &
            gap_thresh, min_length, min_snr, noise_perc, window_size, status, &
            wavelet, omega0) result(peaks)
        complex(dp), intent(in) :: x(:) !! complex input signal
        real(dp), intent(in) :: widths(:) !! positive complex-wavelet scales
        real(dp), intent(in), optional :: max_distances(:) !! maximum ridge shifts
        real(dp), intent(in), optional :: gap_thresh !! maximum missing rows
        integer, intent(in), optional :: min_length !! minimum ridge length
        real(dp), intent(in), optional :: min_snr !! minimum modulus/noise ratio
        real(dp), intent(in), optional :: noise_perc !! local modulus percentile
        integer, intent(in), optional :: window_size !! noise neighborhood length
        integer, intent(out), optional :: status !! success or invalid input
        procedure(cwt_complex_wavelet_callback), optional :: wavelet !! custom complex kernel
        real(dp), intent(in), optional :: omega0 !! default Morlet angular frequency
        integer, allocatable :: peaks(:) !! sorted one-based peak locations

        call complex_peak_kernel(x, widths, peaks, max_distances, gap_thresh, &
            min_length, min_snr, noise_perc, window_size, status, wavelet, omega0)
    end function find_peaks_cwt_complex_complex

    subroutine complex_peak_kernel(x, widths, peaks, max_distances, gap_thresh, &
            min_length, min_snr, noise_perc, window_size, status, wavelet, omega0)
        complex(dp), intent(in) :: x(:) !! real or complex input in complex domain
        real(dp), intent(in) :: widths(:) !! positive wavelet widths
        integer, allocatable, intent(out) :: peaks(:) !! sorted one-based peak locations
        real(dp), intent(in), optional :: max_distances(:) !! maximum ridge shifts
        real(dp), intent(in), optional :: gap_thresh !! maximum missing rows
        integer, intent(in), optional :: min_length !! minimum ridge length
        real(dp), intent(in), optional :: min_snr !! minimum modulus/noise ratio
        real(dp), intent(in), optional :: noise_perc !! local modulus percentile
        integer, intent(in), optional :: window_size !! noise neighborhood length
        integer, intent(out), optional :: status !! success or invalid input
        procedure(cwt_complex_wavelet_callback), optional :: wavelet !! custom complex kernel
        real(dp), intent(in), optional :: omega0 !! Morlet angular frequency
        complex(dp), allocatable :: coefficients(:,:)
        integer :: code

        if (present(status)) status = cwt_peak_invalid_input
        allocate(peaks(0))
        call complex_wavelet_kernel(x, widths, coefficients, wavelet=wavelet, &
            omega0=omega0, status=code)
        if (code /= cwt_peak_success) return
        peaks = find_peaks_cwt_magnitude(coefficients, widths, max_distances, &
            gap_thresh, min_length, min_snr, noise_perc, window_size, status)
    end subroutine complex_peak_kernel

    function ridge_peaks_real(cwt, widths, max_distances, gap_thresh, min_length, &
            min_snr, noise_perc, window_size, status) result(peaks)
        real(dp), intent(in) :: cwt(:,:) !! real wavelet response, scale by sample
        real(dp), intent(in) :: widths(:) !! wavelet widths, positive
        real(dp), intent(in), optional :: max_distances(:) !! maximum per-row ridge displacement
        real(dp), intent(in), optional :: gap_thresh !! maximum missing-scale rows
        integer, intent(in), optional :: min_length !! minimum ridge length
        real(dp), intent(in), optional :: min_snr !! minimum peak/noise ratio
        real(dp), intent(in), optional :: noise_perc !! noise percentile from 0 to 100
        integer, intent(in), optional :: window_size !! noise-window length
        integer, intent(out), optional :: status !! zero on success
        integer, allocatable :: peaks(:) !! one-based peak indices, sorted
        real(dp), allocatable :: distances(:), noises(:)
        integer, allocatable :: firstrow(:), firstcol(:), lastcol(:), gaps(:), lengths(:)
        integer, allocatable :: prev_ids(:), prev_cols(:), temp(:)
        logical, allocatable :: maxima(:,:), active(:)
        real(dp) :: snr_limit, percentile_p, gap_limit, bestdist, delta
        integer :: nw, n, i, j, k, row, startrow, nlines, nprev, closest
        integer :: minlen, win, halfwin, odd, lo, hi, nfound, cap, value

        if (present(status)) status = cwt_peak_invalid_input
        allocate(peaks(0))
        n = size(cwt,2)
        nw = size(widths)
        if (n == 0 .or. nw == 0 .or. size(cwt,1) /= nw) return
        if (any(.not. ieee_is_finite(cwt)) .or. any(.not. ieee_is_finite(widths))) return
        if (any(widths <= 0.0_dp)) return
        allocate(distances(nw))
        distances = widths/4.0_dp
        if (present(max_distances)) then
            if (size(max_distances) < nw) return
            distances = max_distances(:nw)
        end if
        if (any(.not. ieee_is_finite(distances))) return
        if (any(distances < 0.0_dp)) return
        gap_limit = widths(1)
        if (present(gap_thresh)) gap_limit = gap_thresh
        if (.not. ieee_is_finite(gap_limit)) return
        if (gap_limit < 0.0_dp) return
        minlen = ceiling(real(nw,dp)/4.0_dp)
        if (present(min_length)) minlen = min_length
        if (minlen < 1) return
        snr_limit = 1.0_dp
        if (present(min_snr)) snr_limit = min_snr
        if (.not. ieee_is_finite(snr_limit)) return
        if (snr_limit < 0.0_dp) return
        percentile_p = 10.0_dp
        if (present(noise_perc)) percentile_p = noise_perc
        if (.not. ieee_is_finite(percentile_p)) return
        if (percentile_p < 0.0_dp .or. percentile_p > 100.0_dp) return
        win = ceiling(real(n,dp)/20.0_dp)
        if (present(window_size)) win = window_size
        if (win < 1) return

        allocate(maxima(nw,n))
        maxima = .false.
        do row = 1, nw
            do j = 2, n-1
                maxima(row,j) = cwt(row,j) > cwt(row,j-1) .and. &
                    cwt(row,j) > cwt(row,j+1)
            end do
        end do
        startrow = 0
        do row = nw, 1, -1
            if (any(maxima(row,:))) then
                startrow = row
                exit
            end if
        end do
        if (startrow == 0) then
            if (present(status)) status = cwt_peak_success
            return
        end if

        cap = count(maxima)
        allocate(firstrow(cap),firstcol(cap),lastcol(cap),gaps(cap),lengths(cap), &
            active(cap),prev_ids(cap),prev_cols(cap),temp(cap))
        active = .false.
        nlines = 0
        do j = 1, n
            if (.not. maxima(startrow,j)) cycle
            nlines = nlines+1
            firstrow(nlines) = startrow
            firstcol(nlines) = j
            lastcol(nlines) = j
            lengths(nlines) = 1
            gaps(nlines) = 0
            active(nlines) = .true.
        end do
        do row = startrow-1,1,-1
            nprev = 0
            do i = 1,nlines
                if (.not. active(i)) cycle
                nprev = nprev+1
                prev_ids(nprev) = i
                prev_cols(nprev) = lastcol(i)
                gaps(i) = gaps(i)+1
            end do
            do j = 1,n
                if (.not. maxima(row,j)) cycle
                closest = 0
                bestdist = huge(bestdist)
                do k = 1,nprev
                    delta = abs(real(j-prev_cols(k),dp))
                    if (delta < bestdist) then
                        bestdist = delta
                        closest = prev_ids(k)
                    end if
                end do
                if (closest > 0 .and. bestdist <= distances(row)) then
                    firstrow(closest) = row
                    firstcol(closest) = j
                    lastcol(closest) = j
                    lengths(closest) = lengths(closest)+1
                    gaps(closest) = 0
                else
                    nlines = nlines+1
                    firstrow(nlines) = row
                    firstcol(nlines) = j
                    lastcol(nlines) = j
                    lengths(nlines) = 1
                    gaps(nlines) = 0
                    active(nlines) = .true.
                end if
            end do
            do i = 1,nlines
                if (active(i) .and. real(gaps(i),dp) > gap_limit) active(i) = .false.
            end do
        end do

        allocate(noises(n))
        halfwin = win/2
        odd = mod(win,2)
        do j = 1,n
            lo = max(1,j-halfwin)
            hi = min(n,j+halfwin+odd-1)
            ! With even windows SciPy's left-shifted half-open interval
            ! may become empty at the first boundary. Use 0 noise then.
            if (hi < lo) then
                noises(j) = 0.0_dp
            else
                noises(j) = percentile(cwt(1,lo:hi),percentile_p)
            end if
        end do
        nfound = 0
        do i = 1,nlines
            if (lengths(i) < minlen) cycle
            j = firstcol(i)
            if (abs(noises(j)) > 0.0_dp) then
                if (abs(cwt(firstrow(i),j)/noises(j)) < snr_limit) cycle
            end if
            nfound = nfound+1
            temp(nfound) = j
        end do
        ! SciPy returns the accepted ridge columns in ascending order.
        do i = 2,nfound
            value = temp(i)
            j = i-1
            do while (j >= 1)
                if (temp(j) <= value) exit
                temp(j+1) = temp(j)
                j = j-1
            end do
            temp(j+1) = value
        end do
        deallocate(peaks)
        allocate(peaks(nfound))
        peaks = temp(:nfound)
        if (present(status)) status = cwt_peak_success
    end function ridge_peaks_real
end module scifort_signal_cwt_peaks
