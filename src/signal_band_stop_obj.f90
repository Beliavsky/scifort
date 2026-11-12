! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Continuous-order objective used when optimizing the passband edges of
! analog band-stop Butterworth, Chebyshev, and elliptic filters. This is the
! SciPy scipy.signal.band_stop_obj numerical calculation, independently
! implemented from the standard low-pass order equations and complete
! elliptic integrals. The passband index is ZERO based, matching SciPy.
module scifort_signal_band_stop_obj
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_elliptic, only : elliptic_k, elliptic_k_complement
    implicit none
    private

    integer, parameter, public :: band_stop_success = 0
    integer, parameter, public :: band_stop_invalid_input = 1
    real(dp), parameter :: db_to_log = log(10.0_dp)/10.0_dp

    public :: band_stop_obj

contains

    pure real(dp) function log_expm1_positive(x) result(y)
        real(dp), intent(in) :: x !! strictly positive argument, in natural-log units
        ! Preserve accuracy for low ripple (< 0.001 dB) and prevent overflow
        ! for large attenuation. The small-x branch is a Bernoulli expansion.
        if (x < 1.0e-4_dp) then
            y = log(x) + log(1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + x/24.0_dp)))
        else if (x <= log(2.0_dp)) then
            y = log(exp(x) - 1.0_dp)
        else
            y = x + log(1.0_dp-exp(-x))
        end if
    end function log_expm1_positive

    function band_stop_obj(wp, ind, passb, stopb, gpass, gstop, filter_type, status) result(order)
        real(dp), intent(in) :: wp !! candidate analog passband edge (rad/s)
        integer, intent(in) :: ind !! varied edge index, 0 or 1 as in SciPy
        real(dp), intent(in) :: passb(2) !! original analog passband edges
        real(dp), intent(in) :: stopb(2) !! fixed analog stopband edges
        real(dp), intent(in) :: gpass !! positive maximum passband loss (dB)
        real(dp), intent(in) :: gstop !! stopband attenuation (dB), >= gpass
        character(*), intent(in) :: filter_type !! 'butter', 'cheby', or 'ellip'
        integer, intent(out), optional :: status !! 0 for valid finite result; 1 otherwise
        real(dp) :: order
        real(dp) :: edges(2), nat, denom(2), numerator, logratio, m0, m1
        real(dp) :: k0, k0prime, k1, k1prime

        order = ieee_value(0.0_dp, ieee_quiet_nan)
        if (present(status)) status = band_stop_invalid_input
        if (ind /= 0 .and. ind /= 1) return
        if (.not. ieee_is_finite(wp) .or. .not. ieee_is_finite(gpass) .or. &
            .not. ieee_is_finite(gstop)) return
        if (.not. all(ieee_is_finite(passb)) .or. .not. all(ieee_is_finite(stopb))) return
        if (wp <= 0.0_dp .or. any(passb <= 0.0_dp) .or. any(stopb <= 0.0_dp)) return
        if (gpass <= 0.0_dp .or. gstop <= 0.0_dp .or. gpass > gstop) return
        edges = passb
        edges(ind+1) = wp
        if (edges(1) >= edges(2)) return
        denom = stopb**2-edges(1)*edges(2)
        if (any(denom == 0.0_dp)) return
        numerator = edges(1)-edges(2)
        nat = minval(abs(stopb*numerator/denom))
        if (.not. ieee_is_finite(nat) .or. nat <= 0.0_dp .or. nat == 1.0_dp) return

        logratio = log_expm1_positive(db_to_log*gstop)- &
                   log_expm1_positive(db_to_log*gpass)
        select case (trim(filter_type))
        case ('butter')
            ! As in SciPy, a physically non-band-stop trial can yield a
            ! negative objective when nat < 1; optimization uses those trials.
            order = logratio/(2.0_dp*log(nat))
        case ('cheby')
            if (nat <= 1.0_dp) return
            order = acosh(exp(0.5_dp*logratio))/acosh(nat)
        case ('ellip')
            if (nat <= 1.0_dp) return
            if (gpass == gstop) then
                order = 0.0_dp
                if (present(status)) status = band_stop_success
                return
            end if
            ! K(m) and K(1-m) computed with separate AGM inputs to avoid
            ! cancellation when one elliptic parameter approaches unity.
            m0 = 1.0_dp/(nat*nat)
            m1 = exp(-logratio)
            if (m0 <= 0.0_dp .or. m0 >= 1.0_dp .or. &
                m1 <= 0.0_dp .or. m1 >= 1.0_dp) return
            k0 = elliptic_k(m0)
            k0prime = elliptic_k_complement(m0)
            k1 = elliptic_k(m1)
            k1prime = elliptic_k_complement(m1)
            order = (k0/k0prime)*(k1prime/k1)
        case default
            return
        end select
        if (.not. ieee_is_finite(order)) then
            order = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (present(status)) status = band_stop_success
    end function band_stop_obj

end module scifort_signal_band_stop_obj
