! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Frequency-swept signals and Gaussian-modulated pulses in scipy.signal.
module scifort_signal_waveforms
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: waveform_status_success = 0
    integer, parameter, public :: waveform_status_invalid_input = 1
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: degrees_to_radians = pi/180.0_dp

    public :: chirp, chirp_complex, sweep_poly, gausspulse, gausspulse_cutoff

contains

    ! Integrated frequency, in radians, for all four SciPy chirp methods.
    function chirp_phase(t, f0, t1, f1, method, vertex_zero, status) result(phase)
        real(dp), intent(in) :: t(:) !! times at which the chirp is evaluated
        real(dp), intent(in) :: f0 !! instantaneous frequency at time zero
        real(dp), intent(in) :: t1 !! time at which the final frequency is specified
        real(dp), intent(in) :: f1 !! instantaneous frequency at t1
        character(len=*), intent(in), optional :: method !! linear, quadratic, logarithmic, or hyperbolic
        logical, intent(in), optional :: vertex_zero !! quadratic vertex at time zero (default true)
        integer, intent(out) :: status !! zero on success, nonzero on invalid input
        real(dp), allocatable :: phase(:)

        real(dp) :: beta, ratio, sing
        logical :: at_zero
        character(len=:), allocatable :: kind

        allocate(phase(size(t)))
        phase = ieee_value(0.0_dp, ieee_quiet_nan)
        status = waveform_status_invalid_input
        if (.not. all(ieee_is_finite(t))) return
        if (.not. ieee_is_finite(f0)) return
        if (.not. ieee_is_finite(t1)) return
        if (.not. ieee_is_finite(f1)) return
        if (t1 == 0.0_dp) return

        kind = 'linear'
        if (present(method)) kind = trim(method)
        at_zero = .true.
        if (present(vertex_zero)) at_zero = vertex_zero

        select case (kind)
        case ('linear', 'lin', 'li')
            beta = (f1-f0)/t1
            phase = 2.0_dp*pi*(f0*t + 0.5_dp*beta*t*t)
        case ('quadratic', 'quad', 'q')
            beta = (f1-f0)/(t1*t1)
            if (at_zero) then
                phase = 2.0_dp*pi*(f0*t + (beta/3.0_dp)*t*t*t)
            else
                phase = 2.0_dp*pi*(f1*t + (beta/3.0_dp)*((t1-t)**3-t1**3))
            end if
        case ('logarithmic', 'log', 'lo')
            if (.not. ((f0 > 0.0_dp .and. f1 > 0.0_dp) .or. &
                       (f0 < 0.0_dp .and. f1 < 0.0_dp))) return
            if (f0 == f1) then
                phase = 2.0_dp*pi*f0*t
            else
                ratio = f1/f0
                if (.not. ieee_is_finite(ratio)) return
                if (ratio <= 0.0_dp) return
                beta = t1/log(ratio)
                phase = 2.0_dp*pi*beta*f0*(exp(log(ratio)*t/t1)-1.0_dp)
            end if
        case ('hyperbolic', 'hyp')
            if (f0 == 0.0_dp .or. f1 == 0.0_dp) return
            if (f0 == f1) then
                phase = 2.0_dp*pi*f0*t
            else
                sing = -f1*t1/(f0-f1)
                if (.not. ieee_is_finite(sing)) return
                if (sing == 0.0_dp) return
                if (any(t == sing)) return
                phase = 2.0_dp*pi*(-sing*f0)*log(abs(1.0_dp-t/sing))
            end if
        case default
            return
        end select

        if (.not. all(ieee_is_finite(phase))) then
            phase = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        status = waveform_status_success
    end function chirp_phase

    function chirp(t, f0, t1, f1, method, phi, vertex_zero, status) result(y)
        real(dp), intent(in) :: t(:) !! sample times
        real(dp), intent(in) :: f0 !! instantaneous frequency at time zero
        real(dp), intent(in) :: t1 !! time at which f1 is specified
        real(dp), intent(in) :: f1 !! instantaneous frequency at t1
        character(len=*), intent(in), optional :: method !! sweep type (default linear)
        real(dp), intent(in), optional :: phi !! phase offset in degrees (default zero)
        logical, intent(in), optional :: vertex_zero !! quadratic vertex at zero (default true)
        integer, intent(out), optional :: status !! zero on success, nonzero on invalid input
        real(dp), allocatable :: y(:)

        real(dp), allocatable :: phase(:)
        real(dp) :: angle
        integer :: state

        allocate(y(size(t)))
        y = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = waveform_status_invalid_input
        angle = 0.0_dp
        if (present(phi)) angle = phi
        if (.not. ieee_is_finite(angle)) return

        phase = chirp_phase(t, f0, t1, f1, method, vertex_zero, state)
        if (state /= waveform_status_success) return
        if (.not. all(ieee_is_finite(phase + angle*degrees_to_radians))) return
        y = cos(phase + angle*degrees_to_radians)
        if (present(status)) status = waveform_status_success
    end function chirp

    function chirp_complex(t, f0, t1, f1, method, phi, vertex_zero, status) result(y)
        real(dp), intent(in) :: t(:) !! sample times
        real(dp), intent(in) :: f0 !! instantaneous frequency at time zero
        real(dp), intent(in) :: t1 !! time at which f1 is specified
        real(dp), intent(in) :: f1 !! instantaneous frequency at t1
        character(len=*), intent(in), optional :: method !! sweep type (default linear)
        real(dp), intent(in), optional :: phi !! phase offset in degrees (default zero)
        logical, intent(in), optional :: vertex_zero !! quadratic vertex at zero (default true)
        integer, intent(out), optional :: status !! zero on success, nonzero on invalid input
        complex(dp), allocatable :: y(:)

        real(dp), allocatable :: phase(:)
        real(dp) :: angle, nan
        integer :: state

        allocate(y(size(t)))
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        y = cmplx(nan, nan, dp)
        if (present(status)) status = waveform_status_invalid_input
        angle = 0.0_dp
        if (present(phi)) angle = phi
        if (.not. ieee_is_finite(angle)) return

        phase = chirp_phase(t, f0, t1, f1, method, vertex_zero, state)
        if (state /= waveform_status_success) return
        phase = phase + angle*degrees_to_radians
        if (.not. all(ieee_is_finite(phase))) return
        y = cmplx(cos(phase), sin(phase), dp)
        if (present(status)) status = waveform_status_success
    end function chirp_complex

    function sweep_poly(t, poly, phi, status) result(y)
        real(dp), intent(in) :: t(:) !! sample times
        real(dp), intent(in) :: poly(:) !! descending-order instantaneous-frequency polynomial coefficients
        real(dp), intent(in), optional :: phi !! phase offset in degrees (default zero)
        integer, intent(out), optional :: status !! zero on success, nonzero on invalid input
        real(dp), allocatable :: y(:)

        real(dp), allocatable :: phase(:)
        real(dp) :: angle
        integer :: k, n

        allocate(y(size(t)))
        y = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = waveform_status_invalid_input
        if (.not. all(ieee_is_finite(t))) return
        if (.not. all(ieee_is_finite(poly))) return
        angle = 0.0_dp
        if (present(phi)) angle = phi
        if (.not. ieee_is_finite(angle)) return

        n = size(poly)
        allocate(phase(size(t)))
        phase = 0.0_dp
        ! Horner evaluation of the integral with zero constant of integration.
        ! poly(1) multiplies t**(n-1); poly(n) is the constant frequency.
        if (n > 0) then
            phase = poly(1)/real(n,dp)
            do k = 2, n
                phase = phase*t + poly(k)/real(n-k+1,dp)
            end do
            phase = phase*t
        end if
        phase = 2.0_dp*pi*phase + angle*degrees_to_radians
        if (.not. all(ieee_is_finite(phase))) return
        y = cos(phase)
        if (present(status)) status = waveform_status_success
    end function sweep_poly

    ! Envelope factor exp(-a*t*t) as defined in scipy.signal.gausspulse.
    function pulse_factor(fc, bw, bwr, status) result(a)
        real(dp), intent(in), optional :: fc !! center frequency (default 1000)
        real(dp), intent(in), optional :: bw !! fractional bandwidth (default 0.5)
        real(dp), intent(in), optional :: bwr !! bandwidth reference in dB (default -6)
        integer, intent(out) :: status !! zero on success, nonzero on invalid input
        real(dp) :: a, frequency, bandwidth, reference

        a = ieee_value(0.0_dp, ieee_quiet_nan)
        status = waveform_status_invalid_input
        frequency = 1000.0_dp
        bandwidth = 0.5_dp
        reference = -6.0_dp
        if (present(fc)) frequency = fc
        if (present(bw)) bandwidth = bw
        if (present(bwr)) reference = bwr
        if (.not. ieee_is_finite(frequency)) return
        if (.not. ieee_is_finite(bandwidth)) return
        if (.not. ieee_is_finite(reference)) return
        if (frequency < 0.0_dp .or. bandwidth <= 0.0_dp .or. reference >= 0.0_dp) return
        a = -5.0_dp*(pi*frequency*bandwidth)**2/(reference*log(10.0_dp))
        if (.not. ieee_is_finite(a)) return
        status = waveform_status_success
    end function pulse_factor

    function gausspulse(t, fc, bw, bwr, quadrature, envelope, status) result(y)
        real(dp), intent(in) :: t(:) !! sample times
        real(dp), intent(in), optional :: fc !! center frequency (default 1000)
        real(dp), intent(in), optional :: bw !! fractional bandwidth (default 0.5)
        real(dp), intent(in), optional :: bwr !! reference level in dB (default -6)
        real(dp), allocatable, intent(out), optional :: quadrature(:) !! sine component if requested
        real(dp), allocatable, intent(out), optional :: envelope(:) !! pulse envelope if requested
        integer, intent(out), optional :: status !! zero on success, nonzero on invalid input
        real(dp), allocatable :: y(:)

        real(dp), allocatable :: env(:), oscillation(:)
        real(dp) :: a, frequency, nan
        integer :: state

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        allocate(y(size(t)))
        y = nan
        if (present(quadrature)) then
            allocate(quadrature(size(t)))
            quadrature = nan
        end if
        if (present(envelope)) then
            allocate(envelope(size(t)))
            envelope = nan
        end if
        if (present(status)) status = waveform_status_invalid_input
        if (.not. all(ieee_is_finite(t))) return
        a = pulse_factor(fc, bw, bwr, state)
        if (state /= waveform_status_success) return
        frequency = 1000.0_dp
        if (present(fc)) frequency = fc
        allocate(env(size(t)), oscillation(size(t)))
        env = exp(-a*t*t)
        oscillation = 2.0_dp*pi*frequency*t
        if (.not. all(ieee_is_finite(env))) return
        if (.not. all(ieee_is_finite(oscillation))) return
        y = env*cos(oscillation)
        if (present(quadrature)) quadrature = env*sin(oscillation)
        if (present(envelope)) envelope = env
        if (present(status)) status = waveform_status_success
    end function gausspulse

    function gausspulse_cutoff(fc, bw, bwr, tpr, status) result(time)
        real(dp), intent(in), optional :: fc !! center frequency (default 1000)
        real(dp), intent(in), optional :: bw !! fractional bandwidth (default 0.5)
        real(dp), intent(in), optional :: bwr !! reference level in dB (default -6)
        real(dp), intent(in), optional :: tpr !! cutoff level in dB (default -60)
        integer, intent(out), optional :: status !! zero on success, nonzero on invalid input
        real(dp) :: time

        real(dp) :: a, level
        integer :: state

        time = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = waveform_status_invalid_input
        level = -60.0_dp
        if (present(tpr)) level = tpr
        if (.not. ieee_is_finite(level)) return
        if (level >= 0.0_dp) return
        a = pulse_factor(fc, bw, bwr, state)
        if (state /= waveform_status_success) return
        ! An unmodulated, zero-frequency pulse has no finite amplitude cutoff.
        if (a <= 0.0_dp) return
        time = sqrt(-level*log(10.0_dp)/(20.0_dp*a))
        if (.not. ieee_is_finite(time)) then
            time = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (present(status)) status = waveform_status_success
    end function gausspulse_cutoff

end module scifort_signal_waveforms
