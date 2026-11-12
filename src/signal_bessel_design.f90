! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Bessel/Thomson analog and digital IIR filter design.
! Combines the phase/delay/magnitude analog prototype with s-plane frequency
! substitutions and, for digital filters, bilinear s-to-z mapping.
module scifort_signal_bessel_design
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_analog_prototypes, only : analog_prototype_result, &
        besselap, prototype_success
    use scifort_signal_frequency_transforms, only : analog_zpk_transform_result, &
        frequency_transform_success, lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk
    use scifort_signal_bilinear, only : bilinear_zpk
    use scifort_signal_zpk_conversion, only : filter_tf_result, filter_sos_result, &
        zpk2tf, zpk2sos, zpk_conversion_success
    use scifort_signal_design, only : butter_result, design_status_success
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    public :: bessel

    interface bessel
        module procedure bessel_scalar, bessel_vector
    end interface bessel

contains

    function bessel_scalar(order, cutoff, btype, analog, norm, fs) result(out)
        integer, intent(in) :: order !! filter prototype order, 1 through 20
        real(dp), intent(in) :: cutoff !! critical frequency, rad/s or fs units
        character(*), intent(in), optional :: btype !! lowpass (default) or highpass
        logical, intent(in), optional :: analog !! design an analog s-plane filter
        character(*), intent(in), optional :: norm !! phase (default), delay, or mag
        real(dp), intent(in), optional :: fs !! sample rate; digital default is 2
        type(butter_result) :: out

        out = bessel_vector(order, [cutoff], btype, analog, norm, fs)
    end function bessel_scalar

    function bessel_vector(order, cutoff, btype, analog, norm, fs) result(out)
        integer, intent(in) :: order !! filter prototype order, 1 through 20
        real(dp), intent(in) :: cutoff(:) !! one or two critical frequencies
        character(*), intent(in), optional :: btype !! lowpass/highpass/bandpass/bandstop
        logical, intent(in), optional :: analog !! true for continuous-time design
        character(*), intent(in), optional :: norm !! phase (default), delay, or mag
        real(dp), intent(in), optional :: fs !! sample rate for digital design, default 2
        type(butter_result) :: out
        type(analog_prototype_result) :: proto
        type(analog_zpk_transform_result) :: transformed, digital
        type(filter_tf_result) :: tf
        type(filter_sos_result) :: sections
        real(dp), allocatable :: omega(:)
        real(dp) :: rate, center, bandwidth
        character(:), allocatable :: mode, normalization
        logical :: is_analog

        is_analog = .false.
        if (present(analog)) is_analog = analog
        if (order < 1 .or. order > 20) return
        if (is_analog .and. present(fs)) return
        if (size(cutoff) < 1 .or. size(cutoff) > 2) return
        if (.not. all(ieee_is_finite(cutoff))) return
        if (any(cutoff <= 0.0_dp)) return

        mode = 'lowpass'
        if (present(btype)) mode = trim(btype)
        select case (mode)
        case ('low', 'lowpass', 'high', 'highpass')
            if (size(cutoff) /= 1) return
        case ('band', 'bandpass', 'stop', 'bandstop')
            if (size(cutoff) /= 2) return
            if (cutoff(1) >= cutoff(2)) return
        case default
            return
        end select

        normalization = 'phase'
        if (present(norm)) normalization = trim(norm)
        if (normalization /= 'phase' .and. normalization /= 'delay' .and. &
            normalization /= 'mag') return

        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (.not. is_analog) then
            if (.not. ieee_is_finite(rate) .or. rate <= 0.0_dp) return
            if (rate >= huge(rate)/4.0_dp) return
            if (any(cutoff >= rate/2.0_dp)) return
            omega = (2.0_dp*rate)*tan(pi*cutoff/rate)
        else
            omega = cutoff
        end if
        if (.not. all(ieee_is_finite(omega))) return
        if (any(omega <= 0.0_dp)) return

        proto = besselap(order, normalization)
        if (proto%status /= prototype_success) return
        select case (mode)
        case ('low', 'lowpass')
            transformed = lp2lp_zpk(proto%zeros,proto%poles,proto%gain,omega(1))
        case ('high', 'highpass')
            transformed = lp2hp_zpk(proto%zeros,proto%poles,proto%gain,omega(1))
        case ('band', 'bandpass', 'stop', 'bandstop')
            center = sqrt(omega(1))*sqrt(omega(2))
            bandwidth = omega(2)-omega(1)
            if (.not. ieee_is_finite(center) .or. .not. ieee_is_finite(bandwidth)) return
            if (center <= 0.0_dp .or. bandwidth <= 0.0_dp) return
            if (mode == 'band' .or. mode == 'bandpass') then
                transformed = lp2bp_zpk(proto%zeros,proto%poles,proto%gain,center,bandwidth)
            else
                transformed = lp2bs_zpk(proto%zeros,proto%poles,proto%gain,center,bandwidth)
            end if
        end select
        if (transformed%status /= frequency_transform_success) return

        if (is_analog) then
            tf = zpk2tf(transformed%zeros,transformed%poles,transformed%gain)
        else
            digital = bilinear_zpk(transformed%zeros,transformed%poles,transformed%gain,rate)
            if (digital%status /= frequency_transform_success) return
            if (any(abs(digital%poles) >= 1.0_dp)) return
            tf = zpk2tf(digital%zeros,digital%poles,digital%gain)
        end if
        if (tf%status /= zpk_conversion_success) return
        if (.not. nearly_real(tf%b) .or. .not. nearly_real(tf%a)) return
        out%b = real(tf%b,dp)
        out%a = real(tf%a,dp)

        if (.not. is_analog) then
            sections = zpk2sos(digital%zeros,digital%poles,digital%gain)
            if (sections%status /= zpk_conversion_success) return
            if (.not. nearly_real(reshape(sections%sos,[size(sections%sos)]))) return
            out%sos = real(sections%sos,dp)
        end if
        out%status = design_status_success
    end function bessel_vector

    pure logical function nearly_real(values) result(ok)
        complex(dp), intent(in) :: values(:) !! expected conjugate-symmetric polynomial data

        ok = all(ieee_is_finite(real(values,dp))) .and. &
            all(ieee_is_finite(aimag(values)))
        if (.not. ok) return
        ok = maxval(abs(aimag(values))) <= &
            1.0e-8_dp*max(1.0_dp,maxval(abs(real(values,dp))))
    end function nearly_real

end module scifort_signal_bessel_design
