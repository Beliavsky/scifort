! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style one-dimensional spectral analysis, independent implementation.
! Algorithmic inspiration: Jason Christopherson, spectrum (MIT),
! https://github.com/jchristopherson/spectrum, spectrum_periodogram.f90.
! See CODE_PROVENANCE.md and THIRD_PARTY_LICENSES.md.
module scifort_signal_spectral
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_windows_extra, only : window_samples, chebwin, taylor, &
        kaiser_bessel_derived, general_cosine
    use scifort_signal_dpss, only : dpss, dpss_result, dpss_status_success
    implicit none
    private

    integer, parameter, public :: spectral_status_success = 0
    integer, parameter, public :: spectral_status_invalid_input = 1
    real(dp), parameter :: pi = acos(-1.0_dp)

    ! An explicit, reusable plan keeps Bluestein's chirp, transformed kernel,
    ! or the mixed-radix roots, and scratch space across repeated FFTs. Each plan has
    ! mutable scratch storage: do not share one plan across concurrent calls.
    type, public :: fft_plan
        private
        integer :: n = 0, m = 0, algorithm = 0
        complex(dp), allocatable :: chirp(:), kernel_fft(:), scratch(:)
        complex(dp), allocatable :: twiddle(:), work(:)
    contains
        procedure :: initialize => fft_plan_initialize
        procedure :: execute => fft_plan_execute
        procedure :: inverse => fft_plan_inverse
        procedure :: clear => fft_plan_clear
        procedure :: length => fft_plan_length
    end type fft_plan

    type, public :: spectral_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: power(:)
        integer :: status = spectral_status_invalid_input
    end type spectral_result

    type, public :: cross_spectral_result
        real(dp), allocatable :: frequencies(:)
        complex(dp), allocatable :: power(:)
        integer :: status = spectral_status_invalid_input
    end type cross_spectral_result

    type, public :: stft_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        complex(dp), allocatable :: values(:,:)
        integer :: status = spectral_status_invalid_input
    end type stft_result

    type, public :: spectrogram_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: times(:)
        real(dp), allocatable :: power(:,:) !! PSD (default), magnitude, angle or unwrapped phase
        complex(dp), allocatable :: values(:,:) !! complex STFT in mode='complex'; otherwise unallocated
        integer :: status = spectral_status_invalid_input
    end type spectrogram_result

    public :: get_window, general_cosine, chebwin, taylor, kaiser_bessel_derived, &
        gauss_spline, periodogram, welch, csd, coherence, stft, spectrogram
    ! Internal computational kernels shared with the complex spectral layer.
    public :: stft_edge_mode_valid, fill_stft_edges
    interface fill_stft_edges
        module procedure fill_stft_edges_real, fill_stft_edges_complex
    end interface fill_stft_edges
    public :: forward_fft, forward_fft_planned, spectral_setup, median_bias, array_median
    public :: spectral_named_window, spectral_unwrap_phase
    public :: spectral_real_detrend_callback

    ! User callback receives each segment before window multiplication and FFT.
    ! It may replace the data in place, as in SciPy's callable detrend option.
    abstract interface
        subroutine spectral_real_detrend_callback(segment)
            import dp
            real(dp), intent(inout) :: segment(:) !! one unwindowed real segment
        end subroutine spectral_real_detrend_callback
    end interface

contains

    pure logical function stft_edge_mode_valid(mode) result(valid)
        character(*), intent(in) :: mode !! zeros, even, odd, constant, or none
        select case (mode)
        case ('zeros', 'even', 'odd', 'constant', 'none')
            valid = .true.
        case default
            valid = .false.
        end select
    end function stft_edge_mode_valid

    pure subroutine fill_stft_edges_real(x, left, mode, signal)
        real(dp), intent(in) :: x(:) !! finite signal samples, at least left+1 long
        integer, intent(in) :: left !! samples to extend on each side, >= 0
        character(*), intent(in) :: mode !! validated boundary convention
        real(dp), intent(inout) :: signal(:) !! zeroed array, at least size(x)+2*left long
        integer :: i, n

        n = size(x)
        signal(left+1:left+n) = x
        do i = 1, left
            select case (mode)
            case ('even')
                signal(left+1-i) = x(1+i)
                signal(left+n+i) = x(n-i)
            case ('odd')
                signal(left+1-i) = 2.0_dp*x(1) - x(1+i)
                signal(left+n+i) = 2.0_dp*x(n) - x(n-i)
            case ('constant')
                signal(left+1-i) = x(1)
                signal(left+n+i) = x(n)
            case ('zeros', 'none')
                ! Unfilled entries stay zero.
            end select
        end do
    end subroutine fill_stft_edges_real

    pure subroutine fill_stft_edges_complex(x, left, mode, signal)
        complex(dp), intent(in) :: x(:) !! finite signal samples, at least left+1 long
        integer, intent(in) :: left !! samples to extend on each side, >= 0
        character(*), intent(in) :: mode !! validated boundary convention
        complex(dp), intent(inout) :: signal(:) !! zeroed array, at least size(x)+2*left long
        integer :: i, n

        n = size(x)
        signal(left+1:left+n) = x
        do i = 1, left
            select case (mode)
            case ('even')
                signal(left+1-i) = x(1+i)
                signal(left+n+i) = x(n-i)
            case ('odd')
                signal(left+1-i) = 2.0_dp*x(1) - x(1+i)
                signal(left+n+i) = 2.0_dp*x(n) - x(n-i)
            case ('constant')
                signal(left+1-i) = x(1)
                signal(left+n+i) = x(n)
            case ('zeros', 'none')
                ! Unfilled entries stay zero.
            end select
        end do
    end subroutine fill_stft_edges_complex

    pure subroutine spectral_unwrap_phase(phase)
        real(dp), intent(inout) :: phase(:,:) !! wrapped phase, frequency-by-time; unwrapped along frequency
        real(dp) :: current, preceding, delta, adjusted, offset, twopi
        integer :: i, j

        twopi = 2.0_dp*pi
        do j = 1, size(phase,2)
            if (size(phase,1) == 0) cycle
            preceding = phase(1,j)
            offset = 0.0_dp
            do i = 2, size(phase,1)
                current = phase(i,j)
                delta = current-preceding
                adjusted = modulo(delta+pi,twopi)-pi
                ! Match numpy.unwrap at the ambiguous +/-pi boundary.
                if (adjusted == -pi .and. delta > 0.0_dp) adjusted = pi
                if (abs(delta) < pi) adjusted = delta
                offset = offset+adjusted-delta
                phase(i,j) = current+offset
                preceding = current
            end do
        end do
    end subroutine spectral_unwrap_phase

    pure elemental real(dp) function gauss_spline(x, n) result(value)
        real(dp), intent(in) :: x !! spline knot position
        integer, intent(in) :: n !! nonnegative order of the B-spline
        real(dp) :: variance, z

        if (n < 0) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        variance = (real(n, dp) + 1.0_dp) / 12.0_dp
        ! Guard before division: huge(x)/sqrt(2*variance) might overflow.
        if (abs(x) > 40.0_dp * sqrt(2.0_dp * variance)) then
            value = 0.0_dp
        else
            z = x / sqrt(2.0_dp * variance)
            value = exp(-z*z) / sqrt(2.0_dp * pi * variance)
        end if
    end function gauss_spline

    function get_window(name, n, fftbins, status, parameter, parameter2, coefficients) result(w)
        character(*), intent(in) :: name !! SciPy window family name
        integer, intent(in) :: n !! output length
        logical, intent(in), optional :: fftbins !! periodic if true; default true
        integer, intent(out), optional :: status !! zero on success
        real(dp), intent(in), optional :: parameter !! optional shape parameter or DPSS NW
        real(dp), intent(in), optional :: parameter2 !! second shape parameter or exponential center
        real(dp), intent(in), optional :: coefficients(:) !! general cosine weights
        real(dp), allocatable :: w(:)
        type(dpss_result) :: taper
        integer :: code
        logical :: periodic

        periodic = .true.
        if (present(fftbins)) periodic = fftbins
        select case (trim(name))
        case ('general_cosine', 'general cosine')
            if (.not. present(coefficients) .or. present(parameter) .or. present(parameter2)) then
                allocate(w(0))
                code = spectral_status_invalid_input
            else
                w = general_cosine(n, coefficients, sym=.not.periodic, status=code)
            end if
        case ('dpss', 'slepian')
            if (.not. present(parameter) .or. present(parameter2) .or. present(coefficients) .or. n < 1) then
                allocate(w(0))
                code = spectral_status_invalid_input
            else
                taper = dpss(n, parameter, sym=.not.periodic, norm='approximate')
                code = spectral_status_invalid_input
                allocate(w(0))
                if (taper%status == dpss_status_success) then
                    w = taper%tapers(:,1)
                    code = spectral_status_success
                end if
            end if
        case default
            if (present(coefficients)) then
                allocate(w(0))
                code = spectral_status_invalid_input
            else
                w = window_samples(name, n, fftbins=fftbins, parameter=parameter, &
                                   parameter2=parameter2, status=code)
            end if
        end select
        if (present(status)) status = code
    end function get_window

    ! Apply the same parameter/family checks for forward and inverse estimators.
    ! This accepts the SciFort spelling of SciPy get_window tuples.
    function spectral_named_window(name, n, status, parameter, parameter2, &
            coefficients, fftbins) result(w)
        character(*), intent(in) :: name !! named window family
        integer, intent(in) :: n !! required number of window samples
        integer, intent(out) :: status !! success or invalid window/shape parameter
        real(dp), intent(in), optional :: parameter !! first scalar shape parameter
        real(dp), intent(in), optional :: parameter2 !! second scalar shape parameter
        real(dp), intent(in), optional :: coefficients(:) !! general cosine coefficients
        logical, intent(in), optional :: fftbins !! periodic if true (default); symmetric if false
        real(dp), allocatable :: w(:)

        status = spectral_status_invalid_input
        allocate(w(0))
        if (n < 1) return
        if (present(parameter)) then
            select case (trim(name))
            case ('tukey', 'tuk', 'gaussian', 'gauss', 'gss', 'kaiser', 'ksr', &
                  'exponential', 'poisson', 'general_gaussian', &
                  'general gauss', 'general_gauss', 'ggs', &
                  'general_hamming', 'general hamming', &
                  'chebwin', 'cheb', 'taylor', 'taylorwin', &
                  'kaiser_bessel_derived', 'kaiser bessel derived', 'kbd', &
                  'dpss', 'slepian')
            case default
                return
            end select
        end if
        if (present(parameter2)) then
            select case (trim(name))
            case ('general_gaussian', 'general gauss', 'general_gauss', &
                  'ggs', 'taylor', 'taylorwin', 'exponential', 'poisson')
            case default
                return
            end select
        end if
        w = get_window(name, n, fftbins=fftbins, status=status, parameter=parameter, &
            parameter2=parameter2, coefficients=coefficients)
        if (status /= spectral_status_success) return
        if (.not. all(ieee_is_finite(w))) status = spectral_status_invalid_input
    end function spectral_named_window

    ! For real signals with return_onesided=.false., the spectral values
    ! retain FFT order: DC, positive frequencies, then negative frequencies.
    ! The Nyquist bin of an even-length FFT belongs to the negative half.
    pure subroutine full_fft_frequencies(n, rate, frequencies)
        integer, intent(in) :: n !! positive Fourier-transform length
        real(dp), intent(in) :: rate !! sample rate, positive
        real(dp), intent(out) :: frequencies(:) !! FFT-order frequencies, size n
        integer :: k, signed_k

        do k = 0, n-1
            signed_k = k
            if (k > (n-1)/2) signed_k = k-n
            frequencies(k+1) = real(signed_k,dp)*rate/real(n,dp)
        end do
    end subroutine full_fft_frequencies

    ! In-place, unnormalized forward DFT with exp(-2*pi*i*j*k/n).
    ! Keep this legacy entry point unchanged for source/binary compatibility.
    subroutine forward_fft(values)
        complex(dp), intent(inout) :: values(:) !! signal in; Fourier coefficients out
        integer :: n

        n = size(values)
        if (n <= 1) return
        if (iand(n, n - 1) == 0) then
            call radix2_fft(values)
        else if (n <= 16) then
            call direct_fft(values)
        else if (mixed_radix_compatible(n)) then
            call mixed_radix_fft(values)
        else
            call bluestein_fft(values)
        end if
    end subroutine forward_fft

    ! Explicit plan-based alternative; no hidden global cache and no changes
    ! to the original forward_fft interface or to its results.
    subroutine forward_fft_planned(values, plan, status)
        complex(dp), intent(inout) :: values(:) !! signal in; Fourier coefficients out
        type(fft_plan), intent(inout) :: plan !! initialized transform plan
        integer, intent(out), optional :: status !! zero on success, nonzero if invalid
        integer :: code

        call plan%execute(values, code)
        if (present(status)) status = code
    end subroutine forward_fft_planned

    ! Initialize transactionally. Invalid sizes or allocation errors leave any
    ! previously prepared plan unchanged. This is safe even for repeated init.
    subroutine fft_plan_initialize(self, n, status)
        use, intrinsic :: iso_fortran_env, only : int64
        class(fft_plan), intent(inout) :: self !! plan to initialize
        integer, intent(in) :: n !! number of complex samples per transform
        integer, intent(out) :: status !! initialization result
        type(fft_plan) :: trial
        integer :: j, ierr
        integer(int64) :: period, square_mod
        real(dp) :: angle

        status = spectral_status_invalid_input
        if (n < 1) return
        trial%n = n
        if (n == 1 .or. iand(n, n - 1) == 0) then
            trial%algorithm = 1
        else if (n <= 16) then
            trial%algorithm = 2
        else if (mixed_radix_compatible(n)) then
            ! Fast Cooley-Tukey path for products of factors 2, 3, 5, 7.
            ! No per-transform allocation or twiddle recalculation is needed.
            trial%algorithm = 4
            allocate(trial%twiddle(n), trial%scratch(n), trial%work(n), stat=ierr)
            if (ierr /= 0) return
            do j = 0, n-1
                angle = -2.0_dp*pi*real(j,dp)/real(n,dp)
                trial%twiddle(j+1) = cmplx(cos(angle),sin(angle),dp)
            end do
        else
            ! Compute the padded convolution length without integer overflow.
            if (n > (huge(n)-1)/2) return
            trial%algorithm = 3
            trial%m = 1
            do while (trial%m < 2*n-1)
                if (trial%m > huge(n)/2) return
                trial%m = 2*trial%m
            end do
            allocate(trial%chirp(n), trial%kernel_fft(trial%m), &
                trial%scratch(trial%m), stat=ierr)
            if (ierr /= 0) return
            trial%kernel_fft = (0.0_dp,0.0_dp)
            period = 2_int64*int(n,int64)
            do j = 0, n-1
                square_mod = modulo(int(j,int64)*int(j,int64),period)
                angle = pi*real(square_mod,dp)/real(n,dp)
                trial%chirp(j+1) = cmplx(cos(angle),-sin(angle),dp)
                trial%kernel_fft(j+1) = conjg(trial%chirp(j+1))
                if (j > 0) trial%kernel_fft(trial%m-j+1) = trial%kernel_fft(j+1)
            end do
            call radix2_fft(trial%kernel_fft)
            trial%scratch = (0.0_dp,0.0_dp)
        end if
        self%n = trial%n
        self%m = trial%m
        self%algorithm = trial%algorithm
        call move_alloc(trial%chirp, self%chirp)
        call move_alloc(trial%kernel_fft, self%kernel_fft)
        call move_alloc(trial%scratch, self%scratch)
        call move_alloc(trial%twiddle, self%twiddle)
        call move_alloc(trial%work, self%work)
        status = spectral_status_success
    end subroutine fft_plan_initialize

    ! This routine mutates the plan's reusable scratch array. A separate plan
    ! is needed per concurrently executing stream/thread.
    subroutine fft_plan_execute(self, values, status)
        class(fft_plan), intent(inout) :: self !! preinitialized transform plan
        complex(dp), intent(inout) :: values(:) !! input samples, overwritten on success
        integer, intent(out) :: status !! invalid plan or size leaves input unchanged

        status = spectral_status_invalid_input
        if (self%algorithm == 0 .or. size(values) /= self%n) return
        select case (self%algorithm)
        case (1)
            call radix2_fft(values)
        case (2)
            call direct_fft(values)
        case (3)
            self%scratch = (0.0_dp,0.0_dp)
            self%scratch(1:self%n) = values*self%chirp
            call radix2_fft(self%scratch)
            self%scratch = conjg(self%scratch*self%kernel_fft)
            call radix2_fft(self%scratch)
            values = conjg(self%scratch(1:self%n))*self%chirp/real(self%m,dp)
        case (4)
            call mixed_radix_recursive(values, self%scratch, self%work, &
                1, 1, self%n, self%n, self%twiddle)
            values = self%scratch
        case default
            return
        end select
        status = spectral_status_success
    end subroutine fft_plan_execute

    ! Inverse FFT using the same plan (normalized by transform length).
    subroutine fft_plan_inverse(self, values, status)
        class(fft_plan), intent(inout) :: self !! reusable Fourier plan
        complex(dp), intent(inout) :: values(:) !! spectrum in, samples out
        integer, intent(out) :: status !! invalid plan/size leaves values unchanged

        status = spectral_status_invalid_input
        if (self%algorithm == 0 .or. size(values) /= self%n) return
        values = conjg(values)
        call self%execute(values, status)
        values = conjg(values)/real(self%n,dp)
    end subroutine fft_plan_inverse

    subroutine fft_plan_clear(self)
        class(fft_plan), intent(inout) :: self !! storage to release
        self%n = 0
        self%m = 0
        self%algorithm = 0
        if (allocated(self%chirp)) deallocate(self%chirp)
        if (allocated(self%kernel_fft)) deallocate(self%kernel_fft)
        if (allocated(self%scratch)) deallocate(self%scratch)
        if (allocated(self%twiddle)) deallocate(self%twiddle)
        if (allocated(self%work)) deallocate(self%work)
    end subroutine fft_plan_clear

    integer function fft_plan_length(self) result(n)
        class(fft_plan), intent(in) :: self !! plan being inspected
        n = self%n
    end function fft_plan_length

    ! A mixed-radix transform is used only when the entire length factors
    ! into small radices. Lengths with other prime factors use Bluestein.
    pure logical function mixed_radix_compatible(n) result(ok)
        integer, intent(in) :: n !! candidate transform length
        integer :: remaining, r
        integer, parameter :: radices(*) = [2,3,5,7]

        ok = .false.
        if (n < 1) return
        remaining = n
        do r = 1, size(radices)
            do while (mod(remaining, radices(r)) == 0)
                remaining = remaining / radices(r)
            end do
        end do
        ok = remaining == 1
    end function mixed_radix_compatible

    ! Compute a mixed-radix FFT with a one-shot plan. Repeated transforms
    ! should use an explicit fft_plan to reuse roots and workspace.
    subroutine mixed_radix_fft(values)
        complex(dp), intent(inout) :: values(:) !! Fourier input/output
        type(fft_plan) :: plan
        integer :: status

        call plan%initialize(size(values),status)
        if (status == spectral_status_success) then
            call plan%execute(values,status)
        else
            ! Do not silently return the original signal if plan allocation
            ! failed. Bluestein remains the unplanned fallback.
            call bluestein_fft(values)
        end if
    end subroutine mixed_radix_fft

    ! Recursive decimation-in-time Cooley-Tukey FFT. Transform the strided
    ! input subsequences first, then combine m radix-sized butterflies.
    ! The input/output work arrays alternate at each level without allocating.
    !
    ! For k = q + t*m, W_n**(r*k) = W_n**(r*q) * W_radix**(r*t).
    ! Apply W_n**(r*q) once for each q and reuse the radix butterfly for all
    ! t. The first factor is looked up in the full-length root table; because
    ! r < radix and q < m, r*q*(nfft/n) < nfft, so no modulo is required.
    ! Radix 3/5/7 use conjugate-pair sums and differences (forward sign).
    recursive subroutine mixed_radix_recursive(input, output, work, start, stride, n, nfft, twiddle)
        complex(dp), intent(in) :: input(:) !! original, unchanged signal
        complex(dp), intent(inout) :: output(:) !! transformed subsequence
        complex(dp), intent(inout) :: work(:) !! disjoint scratch subsequence
        integer, intent(in) :: start !! one-based input offset
        integer, intent(in) :: stride !! input stride for this level
        integer, intent(in) :: n !! current transform length
        integer, intent(in) :: nfft !! original transform length
        complex(dp), intent(in) :: twiddle(:) !! original-length FFT roots
        real(dp), parameter :: sqrt3half = sqrt(3.0_dp)*0.5_dp
        real(dp), parameter :: c51 = cos(2.0_dp*pi/5.0_dp)
        real(dp), parameter :: c52 = cos(4.0_dp*pi/5.0_dp)
        real(dp), parameter :: s51 = sin(2.0_dp*pi/5.0_dp)
        real(dp), parameter :: s52 = sin(4.0_dp*pi/5.0_dp)
        real(dp), parameter :: c71 = cos(2.0_dp*pi/7.0_dp)
        real(dp), parameter :: c72 = cos(4.0_dp*pi/7.0_dp)
        real(dp), parameter :: c73 = cos(6.0_dp*pi/7.0_dp)
        real(dp), parameter :: s71 = sin(2.0_dp*pi/7.0_dp)
        real(dp), parameter :: s72 = sin(4.0_dp*pi/7.0_dp)
        real(dp), parameter :: s73 = sin(6.0_dp*pi/7.0_dp)
        real(dp), parameter :: cos7(3,3) = reshape([ &
            c71, c72, c73, c72, c73, c71, c73, c71, c72], [3,3])
        real(dp), parameter :: sin7(3,3) = reshape([ &
            s71, s72, s73, s72,-s73,-s71, s73,-s71, s72], [3,3])
        complex(dp), parameter :: imaginary_unit = (0.0_dp,1.0_dp)
        complex(dp) :: a(7), even_part, odd_part, pair_sum(3), pair_diff(3)
        integer :: radix, m, r, q, t, substart, subend, root_step, phase_step

        if (n == 1) then
            output(1) = input(start)
            return
        end if
        if (mod(n,2) == 0) then
            radix = 2
        else if (mod(n,3) == 0) then
            radix = 3
        else if (mod(n,5) == 0) then
            radix = 5
        else
            radix = 7
        end if
        m = n / radix
        do r = 0, radix-1
            substart = r*m+1
            subend = (r+1)*m
            call mixed_radix_recursive(input, work(substart:subend), &
                output(substart:subend), start+r*stride, stride*radix, m, nfft, twiddle)
        end do
        root_step = nfft / n
        do q = 0, m-1
            a(1) = work(q+1)
            phase_step = q*root_step
            do r = 1, radix-1
                a(r+1) = work(r*m+q+1)*twiddle(r*phase_step+1)
            end do
            select case (radix)
            case (2)
                output(q+1) = a(1)+a(2)
                output(q+m+1) = a(1)-a(2)
            case (3)
                pair_sum(1) = a(2)+a(3)
                pair_diff(1) = a(2)-a(3)
                even_part = a(1)-0.5_dp*pair_sum(1)
                odd_part = imaginary_unit*(sqrt3half*pair_diff(1))
                output(q+1) = a(1)+pair_sum(1)
                output(q+m+1) = even_part-odd_part
                output(q+2*m+1) = even_part+odd_part
            case (5)
                pair_sum(1) = a(2)+a(5)
                pair_diff(1) = a(2)-a(5)
                pair_sum(2) = a(3)+a(4)
                pair_diff(2) = a(3)-a(4)
                output(q+1) = a(1)+pair_sum(1)+pair_sum(2)
                even_part = a(1)+c51*pair_sum(1)+c52*pair_sum(2)
                odd_part = imaginary_unit*(s51*pair_diff(1)+s52*pair_diff(2))
                output(q+m+1) = even_part-odd_part
                output(q+4*m+1) = even_part+odd_part
                even_part = a(1)+c52*pair_sum(1)+c51*pair_sum(2)
                odd_part = imaginary_unit*(s52*pair_diff(1)-s51*pair_diff(2))
                output(q+2*m+1) = even_part-odd_part
                output(q+3*m+1) = even_part+odd_part
            case (7)
                do r = 1, 3
                    pair_sum(r) = a(r+1)+a(8-r)
                    pair_diff(r) = a(r+1)-a(8-r)
                end do
                output(q+1) = a(1)+pair_sum(1)+pair_sum(2)+pair_sum(3)
                do t = 1, 3
                    even_part = a(1)
                    odd_part = (0.0_dp,0.0_dp)
                    do r = 1, 3
                        even_part = even_part + cos7(r,t)*pair_sum(r)
                        odd_part = odd_part + sin7(r,t)*pair_diff(r)
                    end do
                    output(q+t*m+1) = even_part-imaginary_unit*odd_part
                    output(q+(7-t)*m+1) = even_part+imaginary_unit*odd_part
                end do
            end select
        end do
    end subroutine mixed_radix_recursive

    ! Small transforms avoid the allocation and padding overhead of Bluestein.
    subroutine direct_fft(values)
        complex(dp), intent(inout) :: values(:) !! short, non-power-of-two signal
        complex(dp) :: temp(size(values))
        real(dp) :: angle
        integer :: j, k, n

        n = size(values)
        do k = 0, n - 1
            temp(k+1) = (0.0_dp, 0.0_dp)
            do j = 0, n - 1
                angle = -2.0_dp*pi*real(k,dp)*real(j,dp)/real(n,dp)
                temp(k+1) = temp(k+1) + values(j+1)*cmplx(cos(angle),sin(angle),dp)
            end do
        end do
        values = temp
    end subroutine direct_fft

    ! Original iterative radix-2 kernel, also used by Bluestein's convolution.
    subroutine radix2_fft(values)
        complex(dp), intent(inout) :: values(:) !! power-of-two-length spectrum
        complex(dp) :: t, w, wm, u
        integer :: n, i, j, k, m, half, stride

        n = size(values)
        if (n <= 1) return
        j = 1
        do i = 1, n - 1
            if (i < j) then
                t = values(i)
                values(i) = values(j)
                values(j) = t
            end if
            k = n / 2
            do while (k >= 1 .and. j > k)
                j = j - k
                k = k / 2
            end do
            j = j + k
        end do
        m = 2
        do while (m <= n)
            half = m / 2
            wm = cmplx(cos(-2.0_dp*pi/real(m,dp)), &
                sin(-2.0_dp*pi/real(m,dp)),dp)
            do stride = 1, n, m
                w = (1.0_dp, 0.0_dp)
                do k = 0, half - 1
                    t = w*values(stride+k+half)
                    u = values(stride+k)
                    values(stride+k) = u + t
                    values(stride+k+half) = u - t
                    w = w*wm
                end do
            end do
            if (m == n) exit
            m = 2*m
        end do
    end subroutine radix2_fft

    ! Chirp convolution: exp(-2*pi*i*j*k/n) =
    ! exp(-pi*i*j**2/n) * exp(-pi*i*k**2/n) * exp(pi*i*(j-k)**2/n).
    ! The convolution is zero-padded to a power of two >= 2*n-1.
    subroutine bluestein_fft(values)
        use, intrinsic :: iso_fortran_env, only : int64
        complex(dp), intent(inout) :: values(:) !! non-power-of-two samples
        complex(dp), allocatable :: a(:), b(:), chirp(:)
        real(dp) :: angle
        integer :: n, m, j
        integer(int64) :: j2mod, period

        n = size(values)
        m = 1
        do while (m < 2*n - 1)
            m = 2*m
        end do
        allocate(a(m), b(m), chirp(n))
        a = (0.0_dp, 0.0_dp)
        b = (0.0_dp, 0.0_dp)
        period = 2_int64*int(n,int64)
        do j = 0, n - 1
            ! Reducing j*j modulo 2*n avoids loss of phase for large j.
            j2mod = modulo(int(j,int64)*int(j,int64),period)
            angle = pi*real(j2mod,dp)/real(n,dp)
            chirp(j+1) = cmplx(cos(angle),-sin(angle),dp)
            a(j+1) = values(j+1)*chirp(j+1)
            b(j+1) = conjg(chirp(j+1))
            if (j > 0) b(m-j+1) = b(j+1)
        end do
        call radix2_fft(a)
        call radix2_fft(b)
        a = conjg(a*b)
        call radix2_fft(a)
        values = conjg(a(:n))*chirp/real(m,dp)
    end subroutine bluestein_fft

    subroutine detrend_segment(seg, method, ok)
        real(dp), intent(inout) :: seg(:) !! unwindowed segment, modified in place
        character(*), intent(in) :: method !! 'constant', 'linear', or 'none'
        logical, intent(out) :: ok !! false for unknown detrending method
        real(dp) :: mid, slope, denom, meanval
        integer :: i, n

        n = size(seg)
        ok = .true.
        select case (trim(method))
        case ('none', 'false')
        case ('constant')
            if (n > 0) seg = seg - sum(seg) / real(n, dp)
        case ('linear')
            if (n <= 0) return
            meanval = sum(seg) / real(n, dp)
            mid = (real(n, dp) + 1.0_dp) * 0.5_dp
            denom = real(n, dp) * (real(n, dp)**2 - 1.0_dp) / 12.0_dp
            slope = 0.0_dp
            if (n > 1) then
                do i = 1, n
                    slope = slope + (real(i, dp) - mid) * (seg(i) - meanval)
                end do
                slope = slope / denom
            end if
            do i = 1, n
                seg(i) = seg(i) - meanval - slope * (real(i, dp) - mid)
            end do
        case default
            ok = .false.
        end select
    end subroutine detrend_segment

    subroutine segment_transform(x, offset, nw, nf, w, trend, ftransform, ok, plan, detrend_callback)
        real(dp), intent(in) :: x(:) !! possibly padded signal
        integer, intent(in) :: offset !! starting one-based sample index
        integer, intent(in) :: nw !! segment size
        integer, intent(in) :: nf !! transform length
        real(dp), intent(in) :: w(:) !! analysis window
        character(*), intent(in) :: trend !! detrending type
        complex(dp), intent(out) :: ftransform(:) !! full FFT, length nf
        logical, intent(out) :: ok !! valid trend and segment
        type(fft_plan), intent(inout), optional :: plan !! reusable FFT plan
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        real(dp) :: local(nw)

        ok = .false.
        if (offset < 1 .or. offset + nw - 1 > size(x)) return
        local = x(offset:offset+nw-1)
        if (present(detrend_callback)) then
            call detrend_callback(local)
            ok = all(ieee_is_finite(local))
        else
            call detrend_segment(local, trend, ok)
        end if
        if (.not. ok) return
        ftransform = (0.0_dp, 0.0_dp)
        ftransform(1:nw) = cmplx(local * w, 0.0_dp, dp)
        if (present(plan)) then
            call forward_fft_planned(ftransform, plan)
        else
            call forward_fft(ftransform)
        end if
    end subroutine segment_transform

    subroutine spectral_setup(n, fs, window, nperseg, noverlap, nfft, scaling, &
        default_window, default_overlap, nw, no, nf, rate, w, factor, ok, &
            window_parameter, window_parameter2, window_coefficients, window_values)
        integer, intent(in) :: n !! available signal length
        real(dp), intent(in), optional :: fs !! samples per second, default 1
        character(*), intent(in), optional :: window !! named periodic window
        integer, intent(in), optional :: nperseg !! segment size, default min(256,n)
        integer, intent(in), optional :: noverlap !! overlap count
        integer, intent(in), optional :: nfft !! transform length
        character(*), intent(in), optional :: scaling !! 'density' or 'spectrum'
        character(*), intent(in) :: default_window !! default window name
        integer, intent(in) :: default_overlap !! negative means nw/2, else nw/divisor
        integer, intent(out) :: nw !! selected segment length
        integer, intent(out) :: no !! selected overlap
        integer, intent(out) :: nf !! selected FFT length
        real(dp), intent(out) :: rate !! effective sample frequency
        real(dp), allocatable, intent(out) :: w(:) !! window weights
        real(dp), intent(out) :: factor !! power-spectrum normalization
        logical, intent(out) :: ok !! all arguments valid
        real(dp), intent(in), optional :: window_parameter !! shape parameter for named window
        real(dp), intent(in), optional :: window_parameter2 !! second shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! general cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! explicit periodic window samples
        integer :: code
        character(:), allocatable :: winname, scalename

        ok = .false.
        rate = 1.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp .or. n < 1) return
        if (present(window_values)) then
            ! In SciPy a numeric window fixes nperseg; it must fit the input.
            if (present(window) .or. present(window_parameter) .or. &
                present(window_parameter2) .or. present(window_coefficients)) return
            nw = size(window_values)
            if (nw < 1 .or. nw > n) return
            if (present(nperseg)) then
                if (nperseg /= nw) return
            end if
            if (.not. all(ieee_is_finite(window_values))) return
            w = window_values
        else
            nw = min(256, n)
            if (present(nperseg)) then
                if (nperseg < 1) return
                nw = min(n, nperseg)
            end if
            winname = default_window
            if (present(window)) winname = trim(window)
            w = spectral_named_window(winname, nw, code, window_parameter, &
                window_parameter2, window_coefficients)
            if (code /= spectral_status_success) return
        end if
        if (.not. all(ieee_is_finite(w))) return
        no = nw / 2
        if (default_overlap > 0) no = nw / default_overlap
        if (present(noverlap)) no = noverlap
        if (no < 0 .or. no >= nw) return
        nf = nw
        if (present(nfft)) nf = nfft
        if (nf < nw) return
        scalename = 'density'
        if (present(scaling)) scalename = trim(scaling)
        select case (scalename)
        case ('density')
            if (sum(w*w) <= tiny(1.0_dp)) return
            factor = 1.0_dp / (rate * sum(w*w))
        case ('spectrum')
            if (abs(sum(w)) <= sqrt(tiny(1.0_dp))) return
            factor = 1.0_dp / sum(w)**2
        case default
            return
        end select
        ok = .true.
    end subroutine spectral_setup

    function csd(x, y, fs, window, nperseg, noverlap, nfft, detrend, scaling, average, &
            window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! first real signal
        real(dp), intent(in) :: y(:) !! second real signal
        real(dp), intent(in), optional :: fs !! sampling rate
        character(*), intent(in), optional :: window !! default 'hann'
        integer, intent(in), optional :: nperseg !! samples per segment
        integer, intent(in), optional :: noverlap !! shared samples
        integer, intent(in), optional :: nfft !! DFT length
        character(*), intent(in), optional :: detrend !! 'constant' (default), 'linear', 'none'
        character(*), intent(in), optional :: scaling !! 'density' (default), 'spectrum'
        character(*), intent(in), optional :: average !! 'mean' (default) or 'median'

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        logical, intent(in), optional :: return_onesided !! default true for real input; false retains all FFT bins
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(cross_spectral_result) :: r
        type(fft_plan) :: plan
        real(dp), allocatable :: w(:), xx(:), yy(:), original(:), imagpart(:)
        complex(dp), allocatable :: fx(:), fy(:), spectra(:,:)
        real(dp) :: rate, factor, correction
        integer :: n, nw, no, nf, segments, hop, k, j, nfout, plan_status
        logical :: ok, onesided
        character(:), allocatable :: trend, avg

        allocate(r%frequencies(0), r%power(0))
        if (size(x) < 1 .or. size(y) < 1) then
            if (size(x) == 0 .and. size(y) == 0) r%status = spectral_status_success
            return
        end if
        if (.not. (all(ieee_is_finite(x)) .and. all(ieee_is_finite(y)))) return
        n = max(size(x), size(y))
        call spectral_setup(n, fs, window, nperseg, noverlap, nfft, scaling, &
            'hann', 0, nw, no, nf, rate, w, factor, ok, &
            window_parameter, window_parameter2, window_coefficients, window_values)
        if (.not. ok) return
        trend = 'constant'
        if (present(detrend)) trend = trim(detrend)
        if (present(detrend_callback) .and. present(detrend)) return
        if (.not. present(detrend_callback)) then
            if (trend /= 'constant' .and. trend /= 'linear' .and. trend /= 'none') return
        end if
        avg = 'mean'
        if (present(average)) avg = trim(average)
        if (avg /= 'mean' .and. avg /= 'median') return
        hop = nw - no
        segments = 1 + (n - nw) / hop
        onesided = .true.
        if (present(return_onesided)) onesided = return_onesided
        nfout = nf
        if (onesided) nfout = nf / 2 + 1
        call plan%initialize(nf, plan_status)
        if (plan_status /= spectral_status_success) return
        allocate(xx(n), yy(n), fx(nf), fy(nf), spectra(nfout,segments))
        xx = 0.0_dp
        yy = 0.0_dp
        xx(1:size(x)) = x
        yy(1:size(y)) = y
        do j = 1, segments
            call segment_transform(xx, 1+(j-1)*hop, nw, nf, w, trend, fx, ok, plan, detrend_callback)
            if (.not. ok) return
            call segment_transform(yy, 1+(j-1)*hop, nw, nf, w, trend, fy, ok, plan, detrend_callback)
            if (.not. ok) return
            spectra(:,j) = conjg(fx(1:nfout)) * fy(1:nfout) * factor
            if (onesided .and. nf > 2) then
                if (mod(nf,2) == 0) then
                    spectra(2:nfout-1,j) = 2.0_dp * spectra(2:nfout-1,j)
                else
                    spectra(2:nfout,j) = 2.0_dp * spectra(2:nfout,j)
                end if
            else if (nf == 1) then
                ! DC is not doubled.
            end if
        end do
        deallocate(r%frequencies, r%power)
        allocate(r%frequencies(nfout), r%power(nfout))
        if (onesided) then
            r%frequencies = [(real(k-1,dp) * rate / real(nf,dp), k=1,nfout)]
        else
            call full_fft_frequencies(nf, rate, r%frequencies)
        end if
        if (avg == 'mean') then
            r%power = sum(spectra, dim=2) / real(segments,dp)
        else
            allocate(original(segments), imagpart(segments))
            correction = median_bias(segments)
            do k = 1, nfout
                original = real(spectra(k,:))
                imagpart = aimag(spectra(k,:))
                r%power(k) = cmplx(array_median(original), array_median(imagpart), dp) / correction
            end do
        end if
        r%status = spectral_status_success
    end function csd

    function periodogram(x, fs, window, nfft, detrend, scaling, &
            window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! real signal
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! default boxcar
        integer, intent(in), optional :: nfft !! DFT size, truncates x when smaller than size(x)
        character(*), intent(in), optional :: detrend !! default constant
        character(*), intent(in), optional :: scaling !! density or spectrum

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        logical, intent(in), optional :: return_onesided !! default true for real input; false retains all FFT bins
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_result) :: r
        type(cross_spectral_result) :: cross
        character(:), allocatable :: win
        integer :: nuse

        allocate(r%frequencies(0), r%power(0))
        if (size(x) == 0) then
            r%status = spectral_status_success
            return
        end if
        nuse = size(x)
        if (present(nfft)) then
            if (nfft < 1) return
            nuse = min(nuse, nfft)
        end if
        if (.not. present(window_values)) win = 'boxcar'
        if (present(window)) win = trim(window)
        cross = csd(x(:nuse), x(:nuse), fs=fs, window=win, nperseg=nuse, &
            noverlap=0, nfft=nfft, detrend=detrend, scaling=scaling, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            return_onesided=return_onesided, detrend_callback=detrend_callback)
        r%frequencies = cross%frequencies
        r%power = real(cross%power)
        r%status = cross%status
    end function periodogram

    function welch(x, fs, window, nperseg, noverlap, nfft, detrend, scaling, average, &
            window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! real signal
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! default hann
        integer, intent(in), optional :: nperseg !! samples per segment
        integer, intent(in), optional :: noverlap !! shared samples
        integer, intent(in), optional :: nfft !! DFT size
        character(*), intent(in), optional :: detrend !! detrending choice
        character(*), intent(in), optional :: scaling !! density or spectrum
        character(*), intent(in), optional :: average !! mean or median

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        logical, intent(in), optional :: return_onesided !! default true for real input; false retains all FFT bins
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_result) :: r
        type(cross_spectral_result) :: cross

        cross = csd(x, x, fs, window, nperseg, noverlap, nfft, detrend, scaling, average, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            return_onesided=return_onesided, detrend_callback=detrend_callback)
        r%frequencies = cross%frequencies
        r%power = real(cross%power)
        r%status = cross%status
    end function welch

    function coherence(x, y, fs, window, nperseg, noverlap, nfft, detrend, &
            window_parameter, window_parameter2, window_coefficients, window_values, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! first real signal
        real(dp), intent(in) :: y(:) !! second real signal
        real(dp), intent(in), optional :: fs !! sample rate
        character(*), intent(in), optional :: window !! default hann
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap
        integer, intent(in), optional :: nfft !! FFT size
        character(*), intent(in), optional :: detrend !! detrending mode

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectral_result) :: r
        type(cross_spectral_result) :: c, pxx, pyy
        integer :: i
        real(dp) :: denom

        allocate(r%frequencies(0), r%power(0))
        c = csd(x, y, fs, window, nperseg, noverlap, nfft, detrend, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
        if (c%status /= spectral_status_success) return
        pxx = csd(x, x, fs, window, nperseg, noverlap, nfft, detrend, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
        pyy = csd(y, y, fs, window, nperseg, noverlap, nfft, detrend, &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            detrend_callback=detrend_callback)
        if (pxx%status /= spectral_status_success .or. pyy%status /= spectral_status_success) return
        r%frequencies = c%frequencies
        deallocate(r%power)
        allocate(r%power(size(c%power)))
        do i = 1, size(r%power)
            denom = real(pxx%power(i)) * real(pyy%power(i))
            r%power(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            if (denom > tiny(1.0_dp)) r%power(i) = abs(c%power(i))**2 / denom
        end do
        r%status = spectral_status_success
    end function coherence

    function stft(x, fs, window, nperseg, noverlap, nfft, detrend, boundary, padded, scaling, &
            window_parameter, window_parameter2, window_coefficients, window_values, return_onesided, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! real signal
        real(dp), intent(in), optional :: fs !! sampling rate
        character(*), intent(in), optional :: window !! default hann
        integer, intent(in), optional :: nperseg !! default min(256,n)
        integer, intent(in), optional :: noverlap !! default half segment
        integer, intent(in), optional :: nfft !! FFT length
        character(*), intent(in), optional :: detrend !! default none
        character(*), intent(in), optional :: boundary !! zeros (default), even, odd, constant, none
        logical, intent(in), optional :: padded !! zero pad trailing data, default true
        character(*), intent(in), optional :: scaling !! spectrum (default) or psd

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        logical, intent(in), optional :: return_onesided !! default true for real input; false retains all FFT bins
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(stft_result) :: r
        type(stft_result) :: temp

        temp = compute_stft(x, fs, window, nperseg, noverlap, nfft, detrend, &
            boundary, padded, scaling, .false., &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            return_onesided=return_onesided, detrend_callback=detrend_callback)
        r = temp
    end function stft

    function spectrogram(x, fs, window, nperseg, noverlap, nfft, detrend, scaling, &
            window_parameter, window_parameter2, window_coefficients, window_values, &
            return_onesided, mode, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! real signal
        real(dp), intent(in), optional :: fs !! sampling rate
        character(*), intent(in), optional :: window !! default Tukey alpha=0.25
        integer, intent(in), optional :: nperseg !! samples per frame
        integer, intent(in), optional :: noverlap !! default nperseg/8
        integer, intent(in), optional :: nfft !! FFT length
        character(*), intent(in), optional :: detrend !! default constant
        character(*), intent(in), optional :: scaling !! density (default) or spectrum

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        logical, intent(in), optional :: return_onesided !! default true for real input; false retains all FFT bins
        character(*), intent(in), optional :: mode !! psd (default), complex, magnitude, angle or phase
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(spectrogram_result) :: r
        type(stft_result) :: temp
        character(:), allocatable :: actual_scale, chosen_mode
        logical :: onesided
        real(dp) :: rate
        integer :: j, last

        chosen_mode = 'psd'
        if (present(mode)) chosen_mode = trim(mode)
        select case (chosen_mode)
        case ('psd', 'complex', 'magnitude', 'angle', 'phase')
        case default
            allocate(r%frequencies(0), r%times(0), r%power(0,0))
            return
        end select
        onesided = .true.
        if (present(return_onesided)) onesided = return_onesided
        rate = 1.0_dp
        if (present(fs)) rate = fs
        actual_scale = 'density'
        if (present(scaling)) actual_scale = trim(scaling)
        temp = compute_stft(x, fs, window=window, nperseg=nperseg, &
            noverlap=noverlap, nfft=nfft, detrend=detrend, boundary='none', &
            padded=.false., scaling=actual_scale, specmode=.true., &
            window_parameter=window_parameter, window_parameter2=window_parameter2, &
            window_coefficients=window_coefficients, window_values=window_values, &
            return_onesided=return_onesided, &
            spectral_amplitude=chosen_mode /= 'psd', detrend_callback=detrend_callback)
        r%frequencies = temp%frequencies
        r%times = temp%times
        if (temp%status == spectral_status_success) then
            if (chosen_mode == 'complex') then
                r%values = temp%values
                allocate(r%power(0,0))
            else
                allocate(r%power(size(temp%values,1),size(temp%values,2)))
                select case (chosen_mode)
                case ('psd')
                    r%power = real(temp%values,dp)
                case ('magnitude')
                    r%power = abs(temp%values)
                case ('angle', 'phase')
                    r%power = atan2(aimag(temp%values),real(temp%values,dp))
                    ! FFTs of real input have strictly real DC and, for even
                    ! FFT sizes, Nyquist bins. SciPy's complex FFT uses -0
                    ! imaginary parts at these bins; its rFFT uses +0.
                    ! Set the principal-angle branch deterministically at
                    ! negative real values before phase unwrapping.
                    do j = 1, size(r%power,2)
                        if (real(temp%values(1,j),dp) < 0.0_dp) then
                            r%power(1,j) = pi
                            if (.not. onesided) r%power(1,j) = -pi
                        else
                            r%power(1,j) = 0.0_dp
                        end if
                    end do
                    last = size(r%power,1)
                    if (.not. onesided) then
                        if (mod(last,2) == 0) then
                            last = last/2+1
                        else
                            last = 0
                        end if
                    end if
                    if (onesided .and. last > 0) then
                        if (abs(temp%frequencies(last)-rate/2.0_dp) > &
                                16.0_dp*epsilon(rate)*max(1.0_dp,rate)) last = 0
                    end if
                    if (last > 1) then
                        do j = 1, size(r%power,2)
                            if (real(temp%values(last,j),dp) < 0.0_dp) then
                                r%power(last,j) = pi
                                if (.not. onesided) r%power(last,j) = -pi
                            else
                                r%power(last,j) = 0.0_dp
                            end if
                        end do
                    end if
                    if (chosen_mode == 'phase') call spectral_unwrap_phase(r%power)
                end select
            end if
        else
            allocate(r%power(0,0))
        end if
        r%status = temp%status
    end function spectrogram

    function compute_stft(x, fs, window, nperseg, noverlap, nfft, detrend, boundary, padded, scaling, specmode, &
            window_parameter, window_parameter2, window_coefficients, window_values, &
            return_onesided, spectral_amplitude, detrend_callback) result(r)
        real(dp), intent(in) :: x(:) !! unpadded real input
        real(dp), intent(in), optional :: fs !! sampling rate
        character(*), intent(in), optional :: window !! named periodic window
        integer, intent(in), optional :: nperseg !! segment length
        integer, intent(in), optional :: noverlap !! overlap length
        integer, intent(in), optional :: nfft !! DFT size
        character(*), intent(in), optional :: detrend !! 'none', 'constant', 'linear'
        character(*), intent(in), optional :: boundary !! zeros, even, odd, constant, or none
        logical, intent(in), optional :: padded !! append enough zeros for complete final frame
        character(*), intent(in), optional :: scaling !! spectrum/psd for STFT; density/spectrum for spectrogram
        logical, intent(in) :: specmode !! true returns real PSD in complex-valued field

        real(dp), intent(in), optional :: window_parameter !! optional Gaussian/Kaiser/Tukey/etc. shape parameter
        real(dp), intent(in), optional :: window_parameter2 !! optional second window shape parameter
        real(dp), intent(in), optional :: window_coefficients(:) !! optional general_cosine coefficients
        real(dp), intent(in), optional :: window_values(:) !! optional caller-supplied window weights
        logical, intent(in), optional :: return_onesided !! default true for real input; false retains all FFT bins
        logical, intent(in), optional :: spectral_amplitude !! use STFT rather than PSD for spectrogram modes
        procedure(spectral_real_detrend_callback), optional :: detrend_callback !! segment callback; excludes detrend=
        type(stft_result) :: r
        type(fft_plan) :: plan
        real(dp), allocatable :: w(:), signal(:), effective_parameter
        complex(dp), allocatable :: ft(:)
        real(dp) :: rate, factor, powerfactor
        integer :: nw, no, nf, nfout, left, nall, hop, nframes, total, j, k, plan_status
        logical :: ok, padmore, onesided, want_amplitude
        character(:), allocatable :: win, trend, edges, scale

        want_amplitude = .false.
        if (present(spectral_amplitude)) want_amplitude = spectral_amplitude
        allocate(r%frequencies(0), r%times(0), r%values(0,0))
        if (size(x) < 1) then
            r%status = spectral_status_success
            return
        end if
        if (.not. all(ieee_is_finite(x))) return
        win = 'hann'
        if (specmode) win = 'tukey'
        if (present(window)) win = trim(window)
        if (present(window_parameter)) effective_parameter = window_parameter
        if (specmode .and. .not. present(window) .and. .not. present(window_parameter) .and. &
            .not. present(window_values) .and. .not. present(window_coefficients)) &
            effective_parameter = 0.25_dp
        scale = 'spectrum'
        if (specmode) scale = 'density'
        if (present(scaling)) scale = trim(scaling)
        ! Spectrogram normalization has power units; STFT has amplitude units.
        if (specmode) then
            call spectral_setup(size(x), fs, window=window, nperseg=nperseg, &
                noverlap=noverlap, nfft=nfft, scaling=scale, default_window=win, &
                default_overlap=8, nw=nw, no=no, nf=nf, rate=rate, &
                w=w, factor=powerfactor, ok=ok, &
                window_parameter=effective_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients, window_values=window_values)
            if (.not. ok) return
            factor = sqrt(powerfactor)
        else
            if (scale == 'spectrum') then
                call spectral_setup(size(x), fs, window=window, nperseg=nperseg, &
                    noverlap=noverlap, nfft=nfft, scaling=scale, default_window=win, &
                    default_overlap=0, nw=nw, no=no, nf=nf, rate=rate, &
                    w=w, factor=powerfactor, ok=ok, &
                window_parameter=effective_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients, window_values=window_values)
            else if (scale == 'psd') then
                call spectral_setup(size(x), fs, window=window, nperseg=nperseg, &
                    noverlap=noverlap, nfft=nfft, scaling='density', default_window=win, &
                    default_overlap=0, nw=nw, no=no, nf=nf, rate=rate, &
                    w=w, factor=powerfactor, ok=ok, &
                window_parameter=effective_parameter, window_parameter2=window_parameter2, &
                window_coefficients=window_coefficients, window_values=window_values)
            else
                return
            end if
            if (.not. ok) return
            factor = sqrt(powerfactor)
        end if
        trend = 'none'
        if (specmode) trend = 'constant'
        if (present(detrend)) trend = trim(detrend)
        if (present(detrend_callback) .and. present(detrend)) return
        if (.not. present(detrend_callback)) then
            if (trend /= 'none' .and. trend /= 'constant' .and. trend /= 'linear') return
        end if
        edges = 'zeros'
        if (specmode) edges = 'none'
        if (present(boundary)) edges = trim(boundary)
        if (.not. stft_edge_mode_valid(edges)) return
        padmore = .true.
        if (specmode) padmore = .false.
        if (present(padded)) padmore = padded
        left = 0
        if (edges /= 'none') left = nw / 2
        nall = size(x) + 2 * left
        hop = nw - no
        if (padmore) then
            nframes = 1 + (nall - nw + hop - 1) / hop
        else
            nframes = 1 + (nall - nw) / hop
        end if
        nframes = max(1, nframes)
        total = nw + (nframes-1) * hop
        call plan%initialize(nf, plan_status)
        if (plan_status /= spectral_status_success) return
        allocate(signal(max(total,nall)), ft(nf))
        signal = 0.0_dp
        call fill_stft_edges(x,left,edges,signal)
        onesided = .true.
        if (present(return_onesided)) onesided = return_onesided
        nfout = nf
        if (onesided) nfout = nf/2 + 1
        deallocate(r%frequencies, r%times, r%values)
        allocate(r%frequencies(nfout), r%times(nframes), r%values(nfout,nframes))
        if (onesided) then
            r%frequencies = [(real(k-1,dp) * rate / real(nf,dp), k=1,nfout)]
        else
            call full_fft_frequencies(nf, rate, r%frequencies)
        end if
        do j = 1, nframes
            r%times(j) = (real((j-1)*hop,dp) + real(nw,dp)/2.0_dp - &
                merge(real(nw,dp)/2.0_dp, 0.0_dp, edges /= 'none')) / rate
            call segment_transform(signal, 1+(j-1)*hop, nw, nf, w, trend, ft, ok, plan, detrend_callback)
            if (.not. ok) return
            if (specmode .and. .not. want_amplitude) then
                r%values(:,j) = cmplx(abs(ft(1:nfout))**2 * powerfactor, 0.0_dp, dp)
                if (onesided .and. nf > 2) then
                    if (mod(nf,2) == 0) then
                        r%values(2:nfout-1,j) = 2.0_dp * r%values(2:nfout-1,j)
                    else
                        r%values(2:nfout,j) = 2.0_dp * r%values(2:nfout,j)
                    end if
                end if
            else
                r%values(:,j) = ft(1:nfout) * factor
            end if
        end do
        r%status = spectral_status_success
    end function compute_stft

    function median_bias(n) result(bias)
        integer, intent(in) :: n !! number of Welch segments
        real(dp) :: bias
        integer :: j
        bias = 1.0_dp
        do j = 1, (n-1)/2
            bias = bias + 1.0_dp / real(2*j+1,dp) - 1.0_dp / real(2*j,dp)
        end do
    end function median_bias

    function array_median(values) result(med)
        real(dp), intent(in) :: values(:) !! segment periodograms
        real(dp) :: med, temp
        real(dp) :: ordered(size(values))
        integer :: i, j, n
        ordered = values
        n = size(ordered)
        do i = 2, n
            temp = ordered(i)
            j = i - 1
            do while (j >= 1)
                if (ordered(j) <= temp) exit
                ordered(j+1) = ordered(j)
                j = j - 1
            end do
            ordered(j+1) = temp
        end do
        med = ordered((n+1)/2)
        if (mod(n,2) == 0) med = (med + ordered(n/2+1)) / 2.0_dp
    end function array_median

end module scifort_signal_spectral
