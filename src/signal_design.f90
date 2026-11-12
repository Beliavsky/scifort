! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! FIR window design, Butterworth/Chebyshev IIR design, and frequency responses.
! Independent implementations of the windowed-sinc, Butterworth prototype,
! frequency transformation, and bilinear transform formulas. See CODE_PROVENANCE.md.
module scifort_signal_design
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_analog_design, only : analog_iir, analog_tf_result, analog_design_success
    use scifort_signal_spectral, only : spectral_named_window, spectral_status_success
    implicit none
    private

    integer, parameter, public :: design_status_success = 0
    integer, parameter, public :: design_status_invalid_input = 1
    integer, parameter, public :: design_status_singular = 2
    real(dp), parameter :: pi = acos(-1.0_dp)

    type, public :: frequency_response
        real(dp), allocatable :: frequencies(:)
        complex(dp), allocatable :: response(:)
        integer :: status = design_status_invalid_input
    end type frequency_response

    type, public :: group_delay_result
        real(dp), allocatable :: frequencies(:)
        real(dp), allocatable :: delay(:)
        integer :: status = design_status_invalid_input
    end type group_delay_result

    type, public :: butter_result
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: a(:)
        real(dp), allocatable :: sos(:,:)
        integer :: status = design_status_invalid_input
    end type butter_result

    type, public :: kaiser_order_result
        integer :: numtaps = 0
        real(dp) :: beta = 0.0_dp
        integer :: status = design_status_invalid_input
    end type kaiser_order_result

    ! SciPy accepts scalar/array cutoffs and logical or string pass_zero.
    ! The optional string overloads retain every existing positional call.
    interface firwin
        module procedure firwin_vector
        module procedure firwin_scalar
        module procedure firwin_vector_logical
        module procedure firwin_scalar_logical
    end interface firwin
    interface freqz
        module procedure freqz_real
        module procedure freqz_complex
    end interface freqz
    interface freqz_zpk
        module procedure freqz_zpk_real_gain
        module procedure freqz_zpk_complex_gain
        module procedure freqz_zpk_real_roots
        module procedure freqz_zpk_real_roots_complex_gain
        module procedure freqz_zpk_real_zeros_real_gain
        module procedure freqz_zpk_real_zeros_complex_gain
        module procedure freqz_zpk_real_poles_real_gain
        module procedure freqz_zpk_real_poles_complex_gain
    end interface freqz_zpk
    interface group_delay
        module procedure group_delay_real
        module procedure group_delay_complex
        module procedure group_delay_sos_real
        module procedure group_delay_sos_complex
        module procedure group_delay_zpk_complex_gain
        module procedure group_delay_zpk_real_gain
        module procedure group_delay_zpk_real_roots
        module procedure group_delay_zpk_real_roots_complex_gain
        module procedure group_delay_zpk_real_zeros_real_gain
        module procedure group_delay_zpk_real_zeros_complex_gain
        module procedure group_delay_zpk_real_poles_real_gain
        module procedure group_delay_zpk_real_poles_complex_gain
    end interface group_delay
    interface group_delay_zpk
        module procedure group_delay_zpk_complex_gain
        module procedure group_delay_zpk_real_gain
        module procedure group_delay_zpk_real_roots
        module procedure group_delay_zpk_real_roots_complex_gain
        module procedure group_delay_zpk_real_zeros_real_gain
        module procedure group_delay_zpk_real_zeros_complex_gain
        module procedure group_delay_zpk_real_poles_real_gain
        module procedure group_delay_zpk_real_poles_complex_gain
    end interface group_delay_zpk
    interface group_delay_sos
        module procedure group_delay_sos_real
        module procedure group_delay_sos_complex
    end interface group_delay_sos
    interface sosfreqz
        module procedure sosfreqz_real
        module procedure sosfreqz_complex
    end interface sosfreqz
    interface butter
        module procedure butter_scalar
        module procedure butter_vector
    end interface butter
    interface cheby1
        module procedure cheby1_scalar
        module procedure cheby1_vector
    end interface cheby1
    public :: firwin, freqz, freqz_zpk, sosfreqz, group_delay, group_delay_sos, group_delay_zpk, butter, cheby1, sos2tf
    public :: kaiser_beta, kaiser_atten, kaiserord

contains

    pure elemental real(dp) function kaiser_beta(attenuation) result(beta)
        real(dp), intent(in) :: attenuation !! required stopband attenuation in dB

        if (attenuation > 50.0_dp) then
            beta = 0.1102_dp * (attenuation - 8.7_dp)
        else if (attenuation > 21.0_dp) then
            beta = 0.5842_dp * (attenuation - 21.0_dp)**0.4_dp + &
                0.07886_dp * (attenuation - 21.0_dp)
        else
            beta = 0.0_dp
        end if
    end function kaiser_beta

    pure elemental real(dp) function kaiser_atten(numtaps, width) result(attenuation)
        integer, intent(in) :: numtaps !! FIR filter length
        real(dp), intent(in) :: width !! transition width normalized to Nyquist

        attenuation = 2.285_dp * (real(numtaps, dp) - 1.0_dp) * pi * width + 7.95_dp
    end function kaiser_atten

    pure function kaiserord(ripple, width) result(design)
        real(dp), intent(in) :: ripple !! desired ripple attenuation in dB, magnitude >= 8
        real(dp), intent(in) :: width !! positive transition width normalized to Nyquist
        type(kaiser_order_result) :: design
        real(dp) :: attenuation, count

        if (.not. ieee_is_finite(ripple) .or. .not. ieee_is_finite(width)) return
        attenuation = abs(ripple)
        if (attenuation < 8.0_dp .or. width <= 0.0_dp) return
        count = (attenuation - 7.95_dp) / (2.285_dp * pi * width) + 1.0_dp
        if (.not. ieee_is_finite(count)) return
        if (count < 1.0_dp .or. count >= real(huge(design%numtaps),dp)) return
        design%numtaps = ceiling(count)
        design%beta = kaiser_beta(attenuation)
        design%status = design_status_success
    end function kaiserord

    pure real(dp) function sinc_pi(x) result(y)
        real(dp), intent(in) :: x !! real argument of sin(pi*x)/(pi*x)
        if (abs(x) < 1.0e-8_dp) then
            y = 1.0_dp - (pi*x)**2 / 6.0_dp
        else
            y = sin(pi*x) / (pi*x)
        end if
    end function sinc_pi

    function firwin_vector(numtaps, cutoff, pass_zero, window, fs, scale, status, &
            window_parameter, width) result(taps)
        integer, intent(in) :: numtaps !! positive FIR length
        real(dp), intent(in) :: cutoff(:) !! strictly increasing frequencies in (0, fs/2)
        character(*), intent(in), optional :: pass_zero !! default lowpass or bandstop
        character(*), intent(in), optional :: window !! symmetric window, default hamming
        real(dp), intent(in), optional :: fs !! sampling rate, default 2 (Nyquist = 1)
        logical, intent(in), optional :: scale !! normalize first passband gain, default true
        integer, intent(out), optional :: status !! success or invalid input
        real(dp), intent(in), optional :: window_parameter !! shape parameter, e.g. Kaiser beta
        real(dp), intent(in), optional :: width !! positive transition width, in units of fs
        real(dp), allocatable :: taps(:)
        real(dp), allocatable :: edges(:), win(:)
        real(dp) :: sample_rate, m, f0, f1, norm_factor, target, omega, beta
        character(:), allocatable :: mode, wname
        integer :: k, band, ncut, st, i
        logical :: first_pass, use_scale, dc_pass, nyquist_pass, current_pass

        allocate(taps(0))
        if (present(status)) status = design_status_invalid_input
        sample_rate = 2.0_dp
        if (present(fs)) sample_rate = fs
        if (.not. ieee_is_finite(sample_rate)) return
        if (sample_rate <= 0.0_dp .or. numtaps < 1) return
        ncut = size(cutoff)
        if (ncut < 1) return
        if (.not. all(ieee_is_finite(cutoff))) return
        if (any(cutoff <= 0.0_dp) .or. any(cutoff >= sample_rate / 2.0_dp)) return
        if (ncut > 1) then
            if (any(cutoff(2:) <= cutoff(:ncut-1))) return
        end if
        mode = 'lowpass'
        if (ncut > 1) mode = 'bandstop'
        if (present(pass_zero)) mode = trim(pass_zero)
        select case (mode)
        case ('lowpass', 'highpass')
            if (ncut /= 1) return
        case ('bandpass', 'bandstop')
            if (ncut < 2) return
        case default
            return
        end select
        dc_pass = mode == 'lowpass' .or. mode == 'bandstop'
        ! A passband ends at Nyquist iff the edge-count parity and the
        ! DC passband selection agree (SciPy's multiband convention).
        nyquist_pass = (mod(ncut,2) == 0) .eqv. dc_pass
        if (nyquist_pass .and. mod(numtaps, 2) == 0) return
        wname = 'hamming'
        if (present(window)) wname = trim(window)
        if (present(width)) then
            ! SciPy's width overrides the explicit window choice.  Convert
            ! physical transition width to the Nyquist-normalized Kaiser beta.
            if (.not. ieee_is_finite(width)) return
            if (width <= 0.0_dp) return
            beta = kaiser_beta(kaiser_atten(numtaps, width/(sample_rate/2.0_dp)))
            if (.not. ieee_is_finite(beta)) return
            win = spectral_named_window('kaiser', numtaps, st, &
                parameter=beta, fftbins=.false.)
        else
            win = spectral_named_window(wname, numtaps, st, &
                parameter=window_parameter, fftbins=.false.)
        end if
        if (st /= spectral_status_success) return
        allocate(edges(ncut+2))
        edges(1) = 0.0_dp
        edges(2:ncut+1) = 2.0_dp * cutoff / sample_rate
        edges(ncut+2) = 1.0_dp
        deallocate(taps)
        allocate(taps(numtaps))
        taps = 0.0_dp
        first_pass = .true.
        target = 0.0_dp
        current_pass = dc_pass
        do band = 1, ncut+1
            f0 = edges(band)
            f1 = edges(band+1)
            if (current_pass) then
                if (first_pass) then
                    if (f0 == 0.0_dp) then
                        target = 0.0_dp
                    else if (f1 == 1.0_dp) then
                        target = 1.0_dp
                    else
                        target = 0.5_dp * (f0 + f1)
                    end if
                    first_pass = .false.
                end if
                do k = 1, numtaps
                    m = real(k-1, dp) - 0.5_dp * real(numtaps-1, dp)
                    taps(k) = taps(k) + f1 * sinc_pi(f1*m) - f0 * sinc_pi(f0*m)
                end do
            end if
            current_pass = .not. current_pass
        end do
        taps = taps * win
        use_scale = .true.
        if (present(scale)) use_scale = scale
        if (use_scale) then
            norm_factor = 0.0_dp
            do i = 1, numtaps
                m = real(i-1, dp) - 0.5_dp * real(numtaps-1, dp)
                omega = pi * target * m
                norm_factor = norm_factor + taps(i) * cos(omega)
            end do
            if (.not. ieee_is_finite(norm_factor)) then
                deallocate(taps)
                allocate(taps(0))
                return
            end if
            if (abs(norm_factor) <= tiny(norm_factor)) then
                deallocate(taps)
                allocate(taps(0))
                return
            end if
            taps = taps / norm_factor
        end if
        if (present(status)) status = design_status_success
    end function firwin_vector

    function firwin_scalar(numtaps, cutoff, pass_zero, window, fs, scale, status, &
            window_parameter, width) result(taps)
        integer, intent(in) :: numtaps !! positive number of FIR coefficients
        real(dp), intent(in) :: cutoff !! scalar cutoff in physical sample-frequency units
        character(*), intent(in), optional :: pass_zero !! lowpass (default) or highpass
        character(*), intent(in), optional :: window !! symmetric FIR window name
        real(dp), intent(in), optional :: fs !! sample frequency; default 2
        logical, intent(in), optional :: scale !! normalize in the first passband
        integer, intent(out), optional :: status !! success or invalid design
        real(dp), intent(in), optional :: window_parameter !! window shape parameter
        real(dp), intent(in), optional :: width !! positive transition width in fs units
        real(dp), allocatable :: taps(:)

        taps = firwin_vector(numtaps, [cutoff], pass_zero, window, fs, scale, &
            status, window_parameter, width)
    end function firwin_scalar

    function firwin_vector_logical(numtaps, cutoff, pass_zero, window, fs, scale, &
            status, window_parameter, width) result(taps)
        integer, intent(in) :: numtaps !! positive number of FIR coefficients
        real(dp), intent(in) :: cutoff(:) !! strictly increasing band edges
        logical, intent(in) :: pass_zero !! true for DC passband, false for DC stopband
        character(*), intent(in), optional :: window !! symmetric FIR window name
        real(dp), intent(in), optional :: fs !! sample frequency; default 2
        logical, intent(in), optional :: scale !! normalize in the first passband
        integer, intent(out), optional :: status !! success or invalid design
        real(dp), intent(in), optional :: window_parameter !! window shape parameter
        real(dp), intent(in), optional :: width !! positive transition width in fs units
        real(dp), allocatable :: taps(:)
        character(:), allocatable :: mode

        if (pass_zero) then
            mode = 'bandstop'
            if (size(cutoff) == 1) mode = 'lowpass'
        else
            mode = 'bandpass'
            if (size(cutoff) == 1) mode = 'highpass'
        end if
        taps = firwin_vector(numtaps, cutoff, mode, window, fs, scale, &
            status, window_parameter, width)
    end function firwin_vector_logical

    function firwin_scalar_logical(numtaps, cutoff, pass_zero, window, fs, scale, &
            status, window_parameter, width) result(taps)
        integer, intent(in) :: numtaps !! positive number of FIR coefficients
        real(dp), intent(in) :: cutoff !! scalar cutoff in physical sample-frequency units
        logical, intent(in) :: pass_zero !! true for lowpass, false for highpass
        character(*), intent(in), optional :: window !! symmetric FIR window name
        real(dp), intent(in), optional :: fs !! sample frequency; default 2
        logical, intent(in), optional :: scale !! normalize in the first passband
        integer, intent(out), optional :: status !! success or invalid design
        real(dp), intent(in), optional :: window_parameter !! window shape parameter
        real(dp), intent(in), optional :: width !! positive transition width in fs units
        real(dp), allocatable :: taps(:)

        taps = firwin_vector_logical(numtaps, [cutoff], pass_zero, window, fs, scale, &
            status, window_parameter, width)
    end function firwin_scalar_logical

    function make_frequencies(worn, whole, fs, include_nyquist, w, ok) result(f)
        integer, intent(in), optional :: worn !! number of equally spaced frequencies, default 512
        logical, intent(in), optional :: whole !! include full unit circle, default false
        real(dp), intent(in), optional :: fs !! sampling rate, default 2*pi
        logical, intent(in), optional :: include_nyquist !! include Nyquist endpoint (half circle)
        real(dp), intent(in), optional :: w(:) !! optional caller-supplied physical frequencies
        logical, intent(out) :: ok !! true on valid frequencies
        real(dp), allocatable :: f(:)
        integer :: n, i
        real(dp) :: rate, endpoint
        logical :: all_circle, include_end

        ok = .false.
        allocate(f(0))
        rate = 2.0_dp * pi
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp) return
        if (present(w)) then
            if (.not. all(ieee_is_finite(w))) return
            f = w
        else
            n = 512
            if (present(worn)) n = worn
            if (n < 1) return
            all_circle = .false.
            if (present(whole)) all_circle = whole
            include_end = .false.
            if (present(include_nyquist)) include_end = include_nyquist
            deallocate(f)
            allocate(f(n))
            endpoint = rate / 2.0_dp
            if (all_circle) endpoint = rate
            do i = 1, n
                if (include_end .and. .not. all_circle .and. n > 1) then
                    f(i) = endpoint * real(i-1,dp) / real(n-1,dp)
                else
                    f(i) = endpoint * real(i-1,dp) / real(n,dp)
                end if
            end do
        end if
        ok = .true.
    end function make_frequencies

    pure complex(dp) function poly_eval(coefs, q) result(y)
        complex(dp), intent(in) :: coefs(:) !! ascending coefficients of delay polynomial
        complex(dp), intent(in) :: q !! exp(-i*omega)
        integer :: i
        y = (0.0_dp, 0.0_dp)
        do i = size(coefs), 1, -1
            y = y*q + coefs(i)
        end do
    end function poly_eval

    function freqz_real(b, a, worn, whole, fs, include_nyquist, w) result(out)
        real(dp), intent(in) :: b(:) !! feedforward coefficients
        real(dp), intent(in), optional :: a(:) !! feedback coefficients, default [1]
        integer, intent(in), optional :: worn !! frequency-bin count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sampling rate, default 2*pi
        logical, intent(in), optional :: include_nyquist !! include endpoint on half-circle
        real(dp), intent(in), optional :: w(:) !! explicit frequencies
        type(frequency_response) :: out
        complex(dp), allocatable :: ac(:)
        if (present(a)) then
            ac = cmplx(a, 0.0_dp, dp)
        else
            ac = [(1.0_dp,0.0_dp)]
        end if
        out = freqz_core(cmplx(b,0.0_dp,dp),ac,worn,whole,fs,include_nyquist,w)
    end function freqz_real

    function freqz_complex(b, a, worn, whole, fs, include_nyquist, w) result(out)
        complex(dp), intent(in) :: b(:) !! complex feedforward coefficients
        complex(dp), intent(in), optional :: a(:) !! feedback coefficients, default [1]
        integer, intent(in), optional :: worn !! frequency-bin count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sampling rate, default 2*pi
        logical, intent(in), optional :: include_nyquist !! include endpoint on half-circle
        real(dp), intent(in), optional :: w(:) !! explicit frequencies
        type(frequency_response) :: out
        complex(dp), allocatable :: ac(:)
        if (present(a)) then
            ac = a
        else
            ac = [(1.0_dp,0.0_dp)]
        end if
        out = freqz_core(b,ac,worn,whole,fs,include_nyquist,w)
    end function freqz_complex

    function freqz_core(b, a, worn, whole, fs, include_nyquist, w) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        integer, intent(in), optional :: worn !! frequency-bin count
        logical, intent(in), optional :: whole !! evaluate complete circle
        real(dp), intent(in), optional :: fs !! sampling rate
        logical, intent(in), optional :: include_nyquist !! include Nyquist endpoint
        real(dp), intent(in), optional :: w(:) !! explicit frequencies
        type(frequency_response) :: out
        complex(dp) :: q, den, numerator
        real(dp) :: rate, angle
        integer :: k
        logical :: ok
        if (size(b) < 1 .or. size(a) < 1) return
        if (.not. all(ieee_is_finite(real(b,dp)))) return
        if (.not. all(ieee_is_finite(aimag(b)))) return
        if (.not. all(ieee_is_finite(real(a,dp)))) return
        if (.not. all(ieee_is_finite(aimag(a)))) return
        if (all(abs(a) == 0.0_dp)) return
        out%frequencies = make_frequencies(worn,whole,fs,include_nyquist,w,ok)
        if (.not. ok) return
        rate = 2.0_dp*pi
        if (present(fs)) rate = fs
        allocate(out%response(size(out%frequencies)))
        out%status = design_status_success
        do k = 1, size(out%response)
            angle = -2.0_dp*pi*out%frequencies(k)/rate
            q = cmplx(cos(angle),sin(angle),dp)
            den = poly_eval(a,q)
            numerator = poly_eval(b,q)
            if (abs(den) <= tiny(1.0_dp)) then
                out%status = design_status_singular
                out%response(k) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
            else
                out%response(k) = numerator/den
            end if
        end do
    end function freqz_core


    function freqz_zpk_real_gain(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! digital filter zeros
        complex(dp), intent(in) :: p(:) !! digital filter poles
        real(dp), intent(in) :: k !! real scalar gain
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! include full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies in units of fs
        type(frequency_response) :: out

        out = freqz_zpk_core(z,p,cmplx(k,0.0_dp,dp),worn,whole,fs,w)
    end function freqz_zpk_real_gain

    function freqz_zpk_complex_gain(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! digital filter zeros
        complex(dp), intent(in) :: p(:) !! digital filter poles
        complex(dp), intent(in) :: k !! complex scalar gain
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! include full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies in units of fs
        type(frequency_response) :: out

        out = freqz_zpk_core(z,p,k,worn,whole,fs,w)
    end function freqz_zpk_complex_gain

    function freqz_zpk_real_roots(z, p, k, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: z(:) !! real digital filter zeros
        real(dp), intent(in) :: p(:) !! real digital filter poles
        real(dp), intent(in) :: k !! real scalar gain
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! include full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies in units of fs
        type(frequency_response) :: out

        out = freqz_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp), &
            cmplx(k,0.0_dp,dp),worn,whole,fs,w)
    end function freqz_zpk_real_roots

    function freqz_zpk_real_roots_complex_gain(z, p, k, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: z(:) !! real digital filter zeros
        real(dp), intent(in) :: p(:) !! real digital filter poles
        complex(dp), intent(in) :: k !! complex scalar gain
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! include full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies in units of fs
        type(frequency_response) :: out

        out = freqz_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp), &
            k,worn,whole,fs,w)
    end function freqz_zpk_real_roots_complex_gain

    function freqz_zpk_real_zeros_real_gain(z, p, k, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: z(:) !! filter zero locations
        complex(dp), intent(in) :: p(:) !! filter pole locations
        real(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! grid count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit physical frequencies
        type(frequency_response) :: out

        out = freqz_zpk_core(cmplx(z,0.0_dp,dp),p,cmplx(k,0.0_dp,dp),worn,whole,fs,w)
    end function freqz_zpk_real_zeros_real_gain

    function freqz_zpk_real_zeros_complex_gain(z, p, k, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: z(:) !! filter zero locations
        complex(dp), intent(in) :: p(:) !! filter pole locations
        complex(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! grid count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit physical frequencies
        type(frequency_response) :: out

        out = freqz_zpk_core(cmplx(z,0.0_dp,dp),p,k,worn,whole,fs,w)
    end function freqz_zpk_real_zeros_complex_gain

    function freqz_zpk_real_poles_real_gain(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! filter zero locations
        real(dp), intent(in) :: p(:) !! filter pole locations
        real(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! grid count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit physical frequencies
        type(frequency_response) :: out

        out = freqz_zpk_core(z,cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp),worn,whole,fs,w)
    end function freqz_zpk_real_poles_real_gain

    function freqz_zpk_real_poles_complex_gain(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! filter zero locations
        real(dp), intent(in) :: p(:) !! filter pole locations
        complex(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! grid count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit physical frequencies
        type(frequency_response) :: out

        out = freqz_zpk_core(z,cmplx(p,0.0_dp,dp),k,worn,whole,fs,w)
    end function freqz_zpk_real_poles_complex_gain

    function freqz_zpk_core(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! zero locations in the z plane
        complex(dp), intent(in) :: p(:) !! pole locations in the z plane
        complex(dp), intent(in) :: k !! scalar response gain
        integer, intent(in), optional :: worn !! frequency count
        logical, intent(in), optional :: whole !! full frequency circle
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequency positions
        type(frequency_response) :: out
        real(dp) :: rate, angle
        complex(dp) :: numerator, denominator, zz
        integer :: i, j
        logical :: ok

        if (.not. all(ieee_is_finite(real(z,dp)))) return
        if (.not. all(ieee_is_finite(aimag(z)))) return
        if (.not. all(ieee_is_finite(real(p,dp)))) return
        if (.not. all(ieee_is_finite(aimag(p)))) return
        if (.not. ieee_is_finite(real(k,dp)) .or. .not. ieee_is_finite(aimag(k))) return
        out%frequencies = make_frequencies(worn,whole,fs,w=w,ok=ok)
        if (.not. ok) return
        rate = 2.0_dp*pi
        if (present(fs)) rate = fs
        allocate(out%response(size(out%frequencies)))
        out%status = design_status_success
        do i = 1, size(out%response)
            angle = 2.0_dp*pi*out%frequencies(i)/rate
            zz = cmplx(cos(angle),sin(angle),dp)
            numerator = k
            denominator = cmplx(1.0_dp,0.0_dp,dp)
            do j = 1, size(z)
                numerator = numerator*(zz-z(j))
            end do
            do j = 1, size(p)
                denominator = denominator*(zz-p(j))
            end do
            if (abs(denominator) <= tiny(1.0_dp)) then
                out%status = design_status_singular
                out%response(i) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan), &
                    ieee_value(0.0_dp,ieee_quiet_nan),dp)
            else
                out%response(i) = numerator/denominator
            end if
        end do
    end function freqz_zpk_core

    function group_delay_zpk_real_gain(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! digital filter zeros
        complex(dp), intent(in) :: p(:) !! digital filter poles
        real(dp), intent(in) :: k !! real scalar gain
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! include full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies in units of fs
        type(group_delay_result) :: out

        out = group_delay_zpk_core(z,p,cmplx(k,0.0_dp,dp),worn,whole,fs,w)
    end function group_delay_zpk_real_gain

    function group_delay_zpk_complex_gain(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! digital filter zeros
        complex(dp), intent(in) :: p(:) !! digital filter poles
        complex(dp), intent(in) :: k !! complex scalar gain
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! include full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies in units of fs
        type(group_delay_result) :: out

        out = group_delay_zpk_core(z,p,k,worn,whole,fs,w)
    end function group_delay_zpk_complex_gain

    function group_delay_zpk_real_roots(z, p, k, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: z(:) !! real digital filter zeros
        real(dp), intent(in) :: p(:) !! real digital filter poles
        real(dp), intent(in) :: k !! real scalar gain
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! include full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies in units of fs
        type(group_delay_result) :: out

        out = group_delay_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp), &
            cmplx(k,0.0_dp,dp),worn,whole,fs,w)
    end function group_delay_zpk_real_roots

    function group_delay_zpk_real_roots_complex_gain(z, p, k, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: z(:) !! real digital filter zeros
        real(dp), intent(in) :: p(:) !! real digital filter poles
        complex(dp), intent(in) :: k !! complex scalar gain
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! include full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies in units of fs
        type(group_delay_result) :: out

        out = group_delay_zpk_core(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp), &
            k,worn,whole,fs,w)
    end function group_delay_zpk_real_roots_complex_gain

    function group_delay_zpk_real_zeros_real_gain(z, p, k, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: z(:) !! filter zero locations
        complex(dp), intent(in) :: p(:) !! filter pole locations
        real(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! grid count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit physical frequencies
        type(group_delay_result) :: out

        out = group_delay_zpk_core(cmplx(z,0.0_dp,dp),p,cmplx(k,0.0_dp,dp),worn,whole,fs,w)
    end function group_delay_zpk_real_zeros_real_gain

    function group_delay_zpk_real_zeros_complex_gain(z, p, k, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: z(:) !! filter zero locations
        complex(dp), intent(in) :: p(:) !! filter pole locations
        complex(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! grid count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit physical frequencies
        type(group_delay_result) :: out

        out = group_delay_zpk_core(cmplx(z,0.0_dp,dp),p,k,worn,whole,fs,w)
    end function group_delay_zpk_real_zeros_complex_gain

    function group_delay_zpk_real_poles_real_gain(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! filter zero locations
        real(dp), intent(in) :: p(:) !! filter pole locations
        real(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! grid count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit physical frequencies
        type(group_delay_result) :: out

        out = group_delay_zpk_core(z,cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp),worn,whole,fs,w)
    end function group_delay_zpk_real_poles_real_gain

    function group_delay_zpk_real_poles_complex_gain(z, p, k, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! filter zero locations
        real(dp), intent(in) :: p(:) !! filter pole locations
        complex(dp), intent(in) :: k !! scalar gain
        integer, intent(in), optional :: worn !! grid count, default 512
        logical, intent(in), optional :: whole !! use full circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit physical frequencies
        type(group_delay_result) :: out

        out = group_delay_zpk_core(z,cmplx(p,0.0_dp,dp),k,worn,whole,fs,w)
    end function group_delay_zpk_real_poles_complex_gain

    ! Stable group delay from roots rather than expanded transfer polynomials.
    ! For H(z)=k*prod(z-zi)/prod(z-pi), the delay in samples is
    !   np-nz - sum Re[zi/(exp(i*w)-zi)] + sum Re[pi/(exp(i*w)-pi)].
    ! This includes the net z-power when the root counts differ. The gain
    ! affects the phase offset, but not its derivative.
    function group_delay_zpk_core(z, p, gain, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: z(:) !! numerator zeros in the z-plane (can be empty)
        complex(dp), intent(in) :: p(:) !! denominator poles in the z-plane (can be empty)
        complex(dp), intent(in) :: gain !! finite nonzero scalar transfer gain
        integer, intent(in), optional :: worn !! grid size, default 512
        logical, intent(in), optional :: whole !! whole frequency circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit grid, in fs units
        type(group_delay_result) :: out
        complex(dp) :: unit, difference
        real(dp) :: rate, angle, nanval, total, correction, term, updated, threshold
        integer :: i, j
        logical :: ok, singular

        if (.not. all(ieee_is_finite(real(z,dp))) .or. &
                .not. all(ieee_is_finite(aimag(z)))) return
        if (.not. all(ieee_is_finite(real(p,dp))) .or. &
                .not. all(ieee_is_finite(aimag(p)))) return
        if (.not. ieee_is_finite(real(gain,dp)) .or. &
                .not. ieee_is_finite(aimag(gain))) return
        if (abs(gain) == 0.0_dp) return
        out%frequencies = make_frequencies(worn, whole, fs, w=w, ok=ok)
        if (.not. ok) then
            deallocate(out%frequencies)
            return
        end if
        rate = 2.0_dp*pi
        if (present(fs)) rate = fs
        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(out%delay(size(out%frequencies)))
        out%delay = nanval
        out%status = design_status_success
        do i = 1, size(out%delay)
            ! Divide by rate first to avoid overflow in 2*pi*frequency.
            angle = (2.0_dp*pi)*(out%frequencies(i)/rate)
            if (.not. ieee_is_finite(angle)) then
                out%status = design_status_singular
                cycle
            end if
            unit = cmplx(cos(angle), sin(angle), dp)
            total = real(size(p)-size(z),dp)
            correction = 0.0_dp
            singular = .false.
            do j = 1, size(z)
                difference = unit-z(j)
                threshold = 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(z(j)))
                if (abs(difference) <= threshold) then
                    singular = .true.
                    exit
                end if
                term = -real(z(j)/difference,dp)
                ! Compensated summation reduces cancellation in long cascades.
                term = term-correction
                updated = total+term
                correction = (updated-total)-term
                total = updated
            end do
            if (singular) then
                out%status = design_status_singular
                cycle
            end if
            do j = 1, size(p)
                difference = unit-p(j)
                threshold = 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(p(j)))
                if (abs(difference) <= threshold) then
                    singular = .true.
                    exit
                end if
                term = real(p(j)/difference,dp)
                term = term-correction
                updated = total+term
                correction = (updated-total)-term
                total = updated
            end do
            if (singular .or. .not. ieee_is_finite(total)) then
                out%status = design_status_singular
            else
                out%delay(i) = total
            end if
        end do
    end function group_delay_zpk_core

    function group_delay_real(b, a, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: b(:) !! numerator coefficients, increasing delay powers
        real(dp), intent(in), optional :: a(:) !! denominator coefficients, default [1]
        integer, intent(in), optional :: worn !! number of bins, default 512
        logical, intent(in), optional :: whole !! use full frequency circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies, overriding worn
        type(group_delay_result) :: out
        complex(dp), allocatable :: ac(:)

        ac = [cmplx(1.0_dp, 0.0_dp, dp)]
        if (present(a)) ac = cmplx(a, 0.0_dp, dp)
        out = group_delay_core(cmplx(b, 0.0_dp, dp), ac, worn, whole, fs, w)
    end function group_delay_real

    function group_delay_complex(b, a, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: b(:) !! numerator coefficients, increasing delay powers
        complex(dp), intent(in), optional :: a(:) !! denominator coefficients, default [1]
        integer, intent(in), optional :: worn !! number of bins, default 512
        logical, intent(in), optional :: whole !! use full frequency circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequencies, overriding worn
        type(group_delay_result) :: out
        complex(dp), allocatable :: ac(:)

        ac = [cmplx(1.0_dp, 0.0_dp, dp)]
        if (present(a)) ac = a
        out = group_delay_core(b, ac, worn, whole, fs, w)
    end function group_delay_complex

    function group_delay_core(b, a, worn, whole, fs, w) result(out)
        ! For H(q)=B(q)/A(q), q=exp(-i*w), the group delay is
        ! real(q*B'(q)/B(q) - q*A'(q)/A(q)).
        ! Horner evaluation avoids finite differences of wrapped phase.
        complex(dp), intent(in) :: b(:) !! numerator coefficients
        complex(dp), intent(in) :: a(:) !! denominator coefficients
        integer, intent(in), optional :: worn !! number of frequency bins
        logical, intent(in), optional :: whole !! use full frequency circle
        real(dp), intent(in), optional :: fs !! sample rate
        real(dp), intent(in), optional :: w(:) !! caller-specified frequencies
        type(group_delay_result) :: out
        complex(dp), allocatable :: bs(:), as(:)
        complex(dp) :: q, bv, av, bm, am
        real(dp) :: rate, angle, nanval, bscale, ascale, btol, atol
        integer :: k, i
        logical :: ok

        if (size(b) == 0 .or. size(a) == 0) return
        if (.not. all(ieee_is_finite(real(b, dp))) .or. &
            .not. all(ieee_is_finite(aimag(b)))) return
        if (.not. all(ieee_is_finite(real(a, dp))) .or. &
            .not. all(ieee_is_finite(aimag(a)))) return
        bscale = maxval(abs(b))
        ascale = maxval(abs(a))
        if (bscale == 0.0_dp .or. ascale == 0.0_dp) return
        bs = b/bscale
        as = a/ascale
        btol = 32.0_dp*epsilon(1.0_dp)*sum(abs(bs))
        atol = 32.0_dp*epsilon(1.0_dp)*sum(abs(as))
        out%frequencies = make_frequencies(worn, whole, fs, w=w, ok=ok)
        if (.not. ok) then
            deallocate(out%frequencies)
            return
        end if
        rate = 2.0_dp*pi
        if (present(fs)) rate = fs
        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(out%delay(size(out%frequencies)))
        out%delay = nanval
        out%status = design_status_success
        do k = 1, size(out%delay)
            angle = -2.0_dp*pi*out%frequencies(k)/rate
            if (.not. ieee_is_finite(angle)) then
                out%status = design_status_singular
                cycle
            end if
            q = cmplx(cos(angle), sin(angle), dp)
            bv = bs(size(bs))
            bm = cmplx(real(size(bs)-1,dp), 0.0_dp, dp)*bv
            do i = size(bs)-1, 1, -1
                bm = bm*q + real(i-1,dp)*bs(i)
                bv = bv*q + bs(i)
            end do
            av = as(size(as))
            am = cmplx(real(size(as)-1,dp), 0.0_dp, dp)*av
            do i = size(as)-1, 1, -1
                am = am*q + real(i-1,dp)*as(i)
                av = av*q + as(i)
            end do
            ! Horner's moment recurrence already includes the q factor.
            if (abs(bv) <= btol .or. abs(av) <= atol) then
                out%status = design_status_singular
                cycle
            end if
            out%delay(k) = real(bm/bv - am/av, dp)
            if (.not. ieee_is_finite(out%delay(k))) then
                out%delay(k) = nanval
                out%status = design_status_singular
            end if
        end do
    end function group_delay_core

    ! Stable group delay of a cascade of second-order sections.  Unlike
    ! group_delay(b,a), this never expands the cascade into high-degree
    ! polynomials, avoiding catastrophic coefficient cancellation.
    ! For q=exp(-i*w), each section contributes
    !   Re[q*B'(q)/B(q) - q*A'(q)/A(q)].
    function group_delay_sos_real(sos, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: sos(:,:) !! normalized real SOS, (nsections,6)
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! use full unit circle, default false
        real(dp), intent(in), optional :: fs !! sampling rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequency grid, in fs units
        type(group_delay_result) :: out

        out = group_delay_sos_core(cmplx(sos, 0.0_dp, dp), worn, whole, fs, w)
    end function group_delay_sos_real

    function group_delay_sos_complex(sos, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: sos(:,:) !! normalized complex SOS, (nsections,6)
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! use full unit circle, default false
        real(dp), intent(in), optional :: fs !! sampling rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequency grid, in fs units
        type(group_delay_result) :: out

        out = group_delay_sos_core(sos, worn, whole, fs, w)
    end function group_delay_sos_complex

    function group_delay_sos_core(sos, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: sos(:,:) !! numerator/denominator triplets, shape (nsec,6)
        integer, intent(in), optional :: worn !! frequency count, default 512
        logical, intent(in), optional :: whole !! full frequency circle, default false
        real(dp), intent(in), optional :: fs !! sample rate, default 2*pi
        real(dp), intent(in), optional :: w(:) !! explicit frequency grid
        type(group_delay_result) :: out
        real(dp), allocatable :: bscale(:), ascale(:)
        complex(dp) :: q, q2, bv, av, bm, am, b0, b1, b2, a0, a1, a2
        real(dp) :: rate, angle, nanval, total, btol, atol
        integer :: j, k, nsec
        logical :: ok, singular

        nsec = size(sos,1)
        if (nsec == 0 .or. size(sos,2) /= 6) return
        if (.not. all(ieee_is_finite(real(sos, dp))) .or. &
                .not. all(ieee_is_finite(aimag(sos)))) return
        if (any(sos(:,4) /= (1.0_dp,0.0_dp))) return
        allocate(bscale(nsec), ascale(nsec))
        do j = 1, nsec
            bscale(j) = maxval(abs(sos(j,1:3)))
            ascale(j) = maxval(abs(sos(j,4:6)))
        end do
        if (any(bscale == 0.0_dp) .or. any(ascale == 0.0_dp)) return
        out%frequencies = make_frequencies(worn, whole, fs, w=w, ok=ok)
        if (.not. ok) then
            deallocate(out%frequencies)
            return
        end if
        rate = 2.0_dp*pi
        if (present(fs)) rate = fs
        nanval = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(out%delay(size(out%frequencies)))
        out%delay = nanval
        out%status = design_status_success
        do k = 1, size(out%delay)
            angle = -2.0_dp*pi*out%frequencies(k)/rate
            if (.not. ieee_is_finite(angle)) then
                out%status = design_status_singular
                cycle
            end if
            q = cmplx(cos(angle), sin(angle), dp)
            q2 = q*q
            total = 0.0_dp
            singular = .false.
            do j = 1, nsec
                ! Scale each second-order section separately.  This handles
                ! cascades with tiny/huge gains without forming H or its
                ! potentially ill-conditioned high-degree polynomials.
                b0 = sos(j,1)/bscale(j)
                b1 = sos(j,2)/bscale(j)
                b2 = sos(j,3)/bscale(j)
                a0 = sos(j,4)/ascale(j)
                a1 = sos(j,5)/ascale(j)
                a2 = sos(j,6)/ascale(j)
                bv = b0 + b1*q + b2*q2
                av = a0 + a1*q + a2*q2
                btol = 32.0_dp*epsilon(1.0_dp)*(abs(b0)+abs(b1)+abs(b2))
                atol = 32.0_dp*epsilon(1.0_dp)*(abs(a0)+abs(a1)+abs(a2))
                if (abs(bv) <= btol .or. abs(av) <= atol) then
                    singular = .true.
                    exit
                end if
                bm = b1*q + 2.0_dp*b2*q2
                am = a1*q + 2.0_dp*a2*q2
                total = total + real(bm/bv - am/av,dp)
            end do
            if (singular .or. .not. ieee_is_finite(total)) then
                out%status = design_status_singular
            else
                out%delay(k) = total
            end if
        end do
    end function group_delay_sos_core

    function sos2tf(sos) result(out)
        ! Multiply the numerator and denominator polynomials of an SOS cascade.
        ! Coefficients are ordered by increasing powers of z**(-1), as in
        ! scipy.signal.sos2tf.  High-order polynomial coefficients can be
        ! ill-conditioned, so filtering should generally retain SOS form.
        real(dp), intent(in) :: sos(:,:) !! sections, shape (nsections, 6)
        type(butter_result) :: out
        real(dp), allocatable :: next_b(:), next_a(:)
        integer :: n, i, j, old_length

        if (size(sos,2) /= 6) return
        if (.not. all(ieee_is_finite(sos))) return
        n = size(sos,1)
        if (n > (huge(n)-1)/2) return
        out%b = [1.0_dp]
        out%a = [1.0_dp]
        do i = 1, n
            old_length = size(out%b)
            allocate(next_b(old_length+2), next_a(old_length+2))
            next_b = 0.0_dp
            next_a = 0.0_dp
            do j = 1, 3
                next_b(j:j+old_length-1) = next_b(j:j+old_length-1) + &
                    sos(i,j)*out%b
                next_a(j:j+old_length-1) = next_a(j:j+old_length-1) + &
                    sos(i,j+3)*out%a
            end do
            call move_alloc(next_b, out%b)
            call move_alloc(next_a, out%a)
            if (.not. all(ieee_is_finite(out%b)) .or. &
                    .not. all(ieee_is_finite(out%a))) then
                out%status = design_status_singular
                return
            end if
        end do
        out%sos = sos
        out%status = design_status_success
    end function sos2tf

    function sosfreqz_real(sos, worn, whole, fs, w) result(out)
        real(dp), intent(in) :: sos(:,:) !! real SOS matrix with shape (nsections,6)
        integer, intent(in), optional :: worn !! frequency-bin count
        logical, intent(in), optional :: whole !! whole circle
        real(dp), intent(in), optional :: fs !! sampling rate
        real(dp), intent(in), optional :: w(:) !! frequencies
        type(frequency_response) :: out
        out = sosfreqz_core(cmplx(sos,0.0_dp,dp),worn,whole,fs,w)
    end function sosfreqz_real

    function sosfreqz_complex(sos, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: sos(:,:) !! complex SOS matrix with shape (nsections,6)
        integer, intent(in), optional :: worn !! frequency-bin count
        logical, intent(in), optional :: whole !! whole circle
        real(dp), intent(in), optional :: fs !! sampling rate
        real(dp), intent(in), optional :: w(:) !! frequencies
        type(frequency_response) :: out
        out = sosfreqz_core(sos,worn,whole,fs,w)
    end function sosfreqz_complex

    function sosfreqz_core(sos, worn, whole, fs, w) result(out)
        complex(dp), intent(in) :: sos(:,:) !! SOS matrix (nsections,6)
        integer, intent(in), optional :: worn !! frequency-bin count
        logical, intent(in), optional :: whole !! whole circle
        real(dp), intent(in), optional :: fs !! sampling rate
        real(dp), intent(in), optional :: w(:) !! frequencies
        type(frequency_response) :: out
        complex(dp) :: q, den, numerator, h
        real(dp) :: rate, angle
        logical :: ok
        integer :: k, j
        if (size(sos,1) < 1 .or. size(sos,2) /= 6) return
        if (.not. all(ieee_is_finite(real(sos,dp)))) return
        if (.not. all(ieee_is_finite(aimag(sos)))) return
        if (any(sos(:,4) /= (1.0_dp,0.0_dp))) return
        out%frequencies = make_frequencies(worn,whole,fs,w=w,ok=ok)
        if (.not. ok) return
        rate = 2.0_dp*pi
        if (present(fs)) rate = fs
        allocate(out%response(size(out%frequencies)))
        out%status = design_status_success
        do k = 1, size(out%response)
            angle = -2.0_dp*pi*out%frequencies(k)/rate
            q = cmplx(cos(angle),sin(angle),dp)
            h = (1.0_dp,0.0_dp)
            do j = 1, size(sos,1)
                den = poly_eval(sos(j,4:6),q)
                numerator = poly_eval(sos(j,1:3),q)
                if (abs(den) <= tiny(1.0_dp)) then
                    out%status = design_status_singular
                    h = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
                    exit
                end if
                h = h*numerator/den
            end do
            out%response(k) = h
        end do
    end function sosfreqz_core

    function butter_scalar(order, cutoff, btype, fs, analog) result(out)
        integer, intent(in) :: order !! positive analog prototype order, up to 20
        real(dp), intent(in) :: cutoff !! digital critical frequency
        character(*), intent(in), optional :: btype !! lowpass (default) or highpass
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        out = butter_vector(order, [cutoff], btype, fs, analog)
    end function butter_scalar

    function butter_vector(order, cutoff, btype, fs, analog) result(out)
        integer, intent(in) :: order !! positive analog prototype order, up to 20
        real(dp), intent(in) :: cutoff(:) !! one or two critical digital frequencies
        character(*), intent(in), optional :: btype !! lowpass, highpass, bandpass, bandstop
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        out = design_iir(order, cutoff, btype, fs, analog=analog)
    end function butter_vector

    function cheby1_scalar(order, rp, cutoff, btype, fs, analog) result(out)
        integer, intent(in) :: order !! Chebyshev type I prototype order, up to 20
        real(dp), intent(in) :: rp !! maximum passband ripple in dB, strictly positive
        real(dp), intent(in) :: cutoff !! critical digital frequency
        character(*), intent(in), optional :: btype !! lowpass (default) or highpass
        real(dp), intent(in), optional :: fs !! sampling frequency, default 2
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        out = design_iir(order, [cutoff], btype, fs, rp, analog)
    end function cheby1_scalar

    function cheby1_vector(order, rp, cutoff, btype, fs, analog) result(out)
        integer, intent(in) :: order !! Chebyshev type I prototype order, up to 20
        real(dp), intent(in) :: rp !! maximum passband ripple in dB, strictly positive
        real(dp), intent(in) :: cutoff(:) !! one or two critical digital frequencies
        character(*), intent(in), optional :: btype !! lowpass, highpass, bandpass, bandstop
        real(dp), intent(in), optional :: fs !! sampling frequency, default 2
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        out = design_iir(order, cutoff, btype, fs, rp, analog)
    end function cheby1_vector

    function design_iir(order, cutoff, btype, fs, ripple, analog) result(out)
        integer, intent(in) :: order !! positive prototype order, at most 20
        real(dp), intent(in) :: cutoff(:) !! one or two critical digital frequencies
        character(*), intent(in), optional :: btype !! lowpass, highpass, bandpass, bandstop
        real(dp), intent(in), optional :: fs !! sampling rate, default 2
        real(dp), intent(in), optional :: ripple !! Chebyshev I passband ripple in dB
        logical, intent(in), optional :: analog !! analog s-domain design, default false
        type(butter_result) :: out
        type(analog_tf_result) :: tf
        complex(dp), allocatable :: analog_poles(:), poles(:), roots(:)
        complex(dp) :: prototype, value, delta, r1, r2, q, numerator, denominator
        real(dp) :: rate, wc1, wc2, bw, w0, ref, factor, mu, eps_ripple, theta
        integer :: k, j, l, n, nsec, section
        logical, allocatable :: used(:)
        character(:), allocatable :: mode
        real(dp), allocatable :: ar(:), br(:), anew(:), bnew(:)

        if (present(analog)) then
            if (analog) then
                ! Analog designs use angular frequency; fs has no meaning.
                if (present(fs)) return
                if (present(ripple)) then
                    tf = analog_iir(order, cutoff, btype, 'cheby1', ripple)
                else
                    tf = analog_iir(order, cutoff, btype, 'butter')
                end if
                if (tf%status /= analog_design_success) return
                out%b = tf%b
                out%a = tf%a
                ! Digital-only second-order sections are left unallocated.
                out%status = design_status_success
                return
            end if
        end if
        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp .or. order < 1 .or. order > 20) return
        mu = 0.0_dp
        if (present(ripple)) then
            if (.not. ieee_is_finite(ripple)) return
            if (ripple <= 0.0_dp .or. ripple > 100.0_dp) return
            eps_ripple = sqrt((10.0_dp**(ripple/10.0_dp)-1.0_dp))
            mu = asinh(1.0_dp/eps_ripple)/real(order,dp)
        end if
        if (.not. all(ieee_is_finite(cutoff))) return
        if (any(cutoff <= 0.0_dp) .or. any(cutoff >= rate/2.0_dp)) return
        mode = 'lowpass'
        if (present(btype)) mode = trim(btype)
        select case (mode)
        case ('lowpass','highpass')
            if (size(cutoff) /= 1) return
        case ('bandpass','bandstop')
            if (size(cutoff) /= 2) return
            if (cutoff(1) >= cutoff(2)) return
        case default
            return
        end select
        wc1 = 2.0_dp * rate * tan(pi * cutoff(1) / rate)
        if (.not. ieee_is_finite(wc1)) return
        wc2 = 0.0_dp
        if (size(cutoff) == 2) then
            wc2 = 2.0_dp * rate * tan(pi * cutoff(2) / rate)
            if (.not. ieee_is_finite(wc2)) return
        end if
        bw = wc2 - wc1
        w0 = sqrt(wc1*wc2)
        n = order
        if (size(cutoff) == 2) n = 2*order
        allocate(analog_poles(n), poles(n), roots(n))
        do k = 0, order-1
            if (present(ripple)) then
                theta = pi*real(2*k+1,dp)/real(2*order,dp)
                prototype = cmplx(-sinh(mu)*sin(theta), cosh(mu)*cos(theta),dp)
            else
                prototype = cmplx(cos(pi*real(2*k+order+1,dp)/real(2*order,dp)), &
                    sin(pi*real(2*k+order+1,dp)/real(2*order,dp)),dp)
            end if
            select case (mode)
            case ('lowpass')
                analog_poles(k+1) = wc1 * prototype
            case ('highpass')
                analog_poles(k+1) = wc1 / prototype
            case ('bandpass')
                value = bw * prototype
                delta = sqrt(value**2 - cmplx(4.0_dp*w0**2,0.0_dp,dp))
                analog_poles(2*k+1) = (value+delta)/2.0_dp
                analog_poles(2*k+2) = (value-delta)/2.0_dp
            case ('bandstop')
                value = bw / prototype
                delta = sqrt(value**2 - cmplx(4.0_dp*w0**2,0.0_dp,dp))
                analog_poles(2*k+1) = (value+delta)/2.0_dp
                analog_poles(2*k+2) = (value-delta)/2.0_dp
            end select
        end do
        poles = (2.0_dp*rate + analog_poles) / (2.0_dp*rate - analog_poles)
        if (.not. all(ieee_is_finite(real(poles,dp)))) return
        if (.not. all(ieee_is_finite(aimag(poles)))) return
        if (any(abs(poles) >= 1.0_dp + 1.0e-10_dp)) return
        select case (mode)
        case ('lowpass')
            roots = cmplx(-1.0_dp,0.0_dp,dp)
            ref = 0.0_dp
        case ('highpass')
            roots = cmplx(1.0_dp,0.0_dp,dp)
            ref = rate/2.0_dp
        case ('bandpass')
            roots(:order) = cmplx(1.0_dp,0.0_dp,dp)
            roots(order+1:) = cmplx(-1.0_dp,0.0_dp,dp)
            ref = rate/pi * atan(w0/(2.0_dp*rate))
        case ('bandstop')
            r1 = (2.0_dp*rate + cmplx(0.0_dp,w0,dp)) / &
                (2.0_dp*rate - cmplx(0.0_dp,w0,dp))
            do k = 1, order
                roots(2*k-1) = r1
                roots(2*k) = conjg(r1)
            end do
            ref = 0.0_dp
        end select
        nsec = (n+1)/2
        allocate(out%sos(nsec,6), used(n))
        out%sos = 0.0_dp
        used = .false.
        do section = 1, nsec
            ! Choose complex poles in conjugate pairs; leave a real pole alone only for odd order.
            k = 0
            if (mod(n,2) == 1 .and. section == nsec) then
                do j = 1, n
                    if (.not. used(j) .and. abs(aimag(poles(j))) <= 1.0e-10_dp) then
                        k = j
                        exit
                    end if
                end do
                if (k == 0) return
                used(k) = .true.
                out%sos(section,4:6) = [1.0_dp, -real(poles(k),dp), 0.0_dp]
                if (mode == 'lowpass') then
                    out%sos(section,1:3) = [1.0_dp,1.0_dp,0.0_dp]
                else
                    out%sos(section,1:3) = [1.0_dp,-1.0_dp,0.0_dp]
                end if
            else
                do j = 1, n
                    if (.not. used(j)) then
                        k = j
                        exit
                    end if
                end do
                if (k == 0) return
                l = 0
                do j = 1, n
                    if (j == k .or. used(j)) cycle
                    if (abs(aimag(poles(k))) < 1.0e-10_dp) then
                        ! A real pole may need pairing with a *different* real pole.
                        if (abs(aimag(poles(j))) > 1.0e-10_dp) cycle
                    else
                        if (abs(poles(j)-conjg(poles(k))) >= 1.0e-8_dp) cycle
                    end if
                    l = j
                    exit
                end do
                if (l == 0) return
                used(k) = .true.
                used(l) = .true.
                out%sos(section,4:6) = [1.0_dp, -real(poles(k)+poles(l),dp), &
                    real(poles(k)*poles(l),dp)]
                select case (mode)
                case ('lowpass')
                    out%sos(section,1:3) = [1.0_dp,2.0_dp,1.0_dp]
                case ('highpass')
                    out%sos(section,1:3) = [1.0_dp,-2.0_dp,1.0_dp]
                case ('bandpass')
                    out%sos(section,1:3) = [1.0_dp,0.0_dp,-1.0_dp]
                case ('bandstop')
                    out%sos(section,1:3) = [1.0_dp, -2.0_dp*real(roots(1),dp), 1.0_dp]
                end select
            end if
            q = cmplx(cos(-2.0_dp*pi*ref/rate),sin(-2.0_dp*pi*ref/rate),dp)
            numerator = poly_eval(cmplx(out%sos(section,1:3),0.0_dp,dp),q)
            denominator = poly_eval(cmplx(out%sos(section,4:6),0.0_dp,dp),q)
            if (abs(numerator) <= tiny(1.0_dp)) return
            factor = abs(denominator/numerator)
            if (.not. ieee_is_finite(factor)) return
            ! Even-order Chebyshev I prototypes have -rp dB gain at DC.
            if (present(ripple)) then
                if (mod(order,2) == 0 .and. section == 1) then
                    factor = factor * 10.0_dp**(-ripple/20.0_dp)
                end if
            end if
            out%sos(section,1:3) = factor * out%sos(section,1:3)
        end do
        allocate(out%b(1),out%a(1))
        out%b = 1.0_dp
        out%a = 1.0_dp
        do section = 1, nsec
            if (mod(n,2) == 1 .and. section == nsec) then
                br = out%sos(section,1:2)
                ar = out%sos(section,4:5)
            else
                br = out%sos(section,1:3)
                ar = out%sos(section,4:6)
            end if
            bnew = conv_real(out%b,br)
            anew = conv_real(out%a,ar)
            call move_alloc(bnew,out%b)
            call move_alloc(anew,out%a)
        end do
        out%status = design_status_success
    end function design_iir

    pure function conv_real(a,b) result(c)
        real(dp), intent(in) :: a(:) !! first polynomial
        real(dp), intent(in) :: b(:) !! second polynomial
        real(dp) :: c(size(a)+size(b)-1)
        integer :: i,j
        c = 0.0_dp
        do j = 1, size(b)
            do i = 1, size(a)
                c(i+j-1) = c(i+j-1)+a(i)*b(j)
            end do
        end do
    end function conv_real
end module scifort_signal_design
