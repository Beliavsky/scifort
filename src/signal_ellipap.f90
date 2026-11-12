! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Expose the existing elliptic lowpass design as a normalized ZPK prototype.
module scifort_signal_ellipap
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_elliptic, only : elliptic_analog_zpk, elliptic_zpk_result
    use scifort_signal_analog_prototypes, only : analog_prototype_result, &
        prototype_success, prototype_invalid_input, prototype_no_convergence, &
        prototype_too_expensive
    implicit none
    private
    public :: ellipap
contains
    function ellipap(n, rp, rs) result(out)
        integer, intent(in) :: n !! prototype order; supported orders 0 to 20
        real(dp), intent(in) :: rp !! positive passband ripple in decibels
        real(dp), intent(in) :: rs !! stopband attenuation in decibels
        type(analog_prototype_result) :: out
        type(elliptic_zpk_result) :: elliptic
        if (n < 0) return
        if (n > 20) then
            out%status = prototype_too_expensive
            return
        end if
        if (.not. ieee_is_finite(rp) .or. .not. ieee_is_finite(rs)) return
        if (rp <= 0.0_dp .or. rs <= rp .or. rs > 100.0_dp) return
        if (n == 0) then
            allocate(out%zeros(0),out%poles(0))
            out%gain = 10.0_dp**(-rp/20.0_dp)
            out%status = prototype_success
            return
        end if
        elliptic = elliptic_analog_zpk(n,rp,rs,[1.0_dp],'lowpass')
        if (elliptic%status /= 0) then
            out%status = prototype_no_convergence
            return
        end if
        out%zeros = elliptic%zeros
        out%poles = elliptic%poles
        out%gain = elliptic%gain
        out%status = prototype_success
    end function ellipap
end module scifort_signal_ellipap
