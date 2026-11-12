! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Analog lowpass frequency substitutions in transfer-function and ZPK form.
! Polynomial transformations use independent rational-polynomial composition,
! without converting to roots. All coefficient vectors use descending powers.
module scifort_signal_frequency_transforms
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: frequency_transform_success = 0
    integer, parameter, public :: frequency_transform_invalid_input = 1
    integer, parameter, public :: frequency_transform_too_expensive = 2
    integer, parameter :: max_polynomial_order = 32

    type, public :: analog_tf_transform_result
        complex(dp), allocatable :: b(:) !! numerator, descending powers of s
        complex(dp), allocatable :: a(:) !! denominator, descending powers of s
        integer :: status = frequency_transform_invalid_input !! zero on success
    end type analog_tf_transform_result

    type, public :: analog_zpk_transform_result
        complex(dp), allocatable :: zeros(:) !! transformed finite zeros
        complex(dp), allocatable :: poles(:) !! transformed finite poles
        complex(dp) :: gain = (0.0_dp, 0.0_dp) !! transformed analog gain
        integer :: status = frequency_transform_invalid_input !! zero on success
    end type analog_zpk_transform_result

    public :: lp2lp, lp2hp, lp2bp, lp2bs
    public :: lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk

    interface lp2lp
        module procedure lp2lp_real, lp2lp_complex
    end interface
    interface lp2hp
        module procedure lp2hp_real, lp2hp_complex
    end interface
    interface lp2bp
        module procedure lp2bp_real, lp2bp_complex
    end interface
    interface lp2bs
        module procedure lp2bs_real, lp2bs_complex
    end interface
    interface lp2lp_zpk
        module procedure lp2lp_zpk_real_gain, lp2lp_zpk_complex_gain
    end interface
    interface lp2hp_zpk
        module procedure lp2hp_zpk_real_gain, lp2hp_zpk_complex_gain
    end interface
    interface lp2bp_zpk
        module procedure lp2bp_zpk_real_gain, lp2bp_zpk_complex_gain
    end interface
    interface lp2bs_zpk
        module procedure lp2bs_zpk_real_gain, lp2bs_zpk_complex_gain
    end interface

contains

    pure logical function finite_complex(z) result(ok)
        complex(dp), intent(in) :: z(:) !! complex values
        ok = all(ieee_is_finite(real(z, dp))) .and. all(ieee_is_finite(aimag(z)))
    end function finite_complex

    pure function poly_product(x, y) result(p)
        complex(dp), intent(in) :: x(:) !! ascending polynomial coefficients
        complex(dp), intent(in) :: y(:) !! ascending polynomial coefficients
        complex(dp) :: p(size(x)+size(y)-1)
        integer :: i, j

        p = (0.0_dp, 0.0_dp)
        do j = 1, size(y)
            do i = 1, size(x)
                p(i+j-1) = p(i+j-1) + x(i)*y(j)
            end do
        end do
    end function poly_product

    pure function poly_add(x, y) result(p)
        complex(dp), intent(in) :: x(:) !! first ascending polynomial
        complex(dp), intent(in) :: y(:) !! second ascending polynomial
        complex(dp) :: p(max(size(x),size(y)))
        p = (0.0_dp, 0.0_dp)
        p(:size(x)) = x
        p(:size(y)) = p(:size(y)) + y
    end function poly_add

    pure function rational_polynomial(c, degree, num, den) result(p)
        complex(dp), intent(in) :: c(:) !! polynomial coefficients, descending
        integer, intent(in) :: degree !! common rational-denominator exponent
        complex(dp), intent(in) :: num(:) !! ascending polynomial for numerator of t(s)
        complex(dp), intent(in) :: den(:) !! ascending polynomial for denominator of t(s)
        complex(dp), allocatable :: p(:), d_power(:)
        integer :: j

        ! Horner's scheme applied to c(t), clearing denominator powers.
        p = [c(1)]
        d_power = [cmplx(1.0_dp, 0.0_dp, dp)]
        do j = 2, size(c)
            d_power = poly_product(d_power, den)
            p = poly_add(poly_product(p, num), c(j)*d_power)
        end do
        ! Align numerator and denominator to the same power of den(s).
        do j = size(c), degree
            p = poly_product(p, den)
        end do
    end function rational_polynomial

    function transform_ba(b, a, which, wo, bw) result(out)
        complex(dp), intent(in) :: b(:) !! lowpass numerator coefficients
        complex(dp), intent(in) :: a(:) !! lowpass denominator coefficients
        integer, intent(in) :: which !! 1 lowpass, 2 highpass, 3 bandpass, 4 bandstop
        real(dp), intent(in), optional :: wo !! target angular center frequency
        real(dp), intent(in), optional :: bw !! target angular bandwidth
        type(analog_tf_transform_result) :: out
        complex(dp), allocatable :: num(:), den(:), bc(:), ac(:), pb(:), pa(:)
        complex(dp) :: scale
        real(dp) :: omega, bandwidth
        integer :: ib, ia, bdeg, adeg, degree

        omega = 1.0_dp
        bandwidth = 1.0_dp
        if (present(wo)) omega = wo
        if (present(bw)) bandwidth = bw
        if (.not. ieee_is_finite(omega) .or. .not. ieee_is_finite(bandwidth)) return
        if (omega <= 0.0_dp .or. bandwidth <= 0.0_dp) return
        ! Band transformations explicitly form wo**2; reject an unrepresentable
        ! center-frequency square before evaluating it under FPE traps.
        if (which >= 3) then
            if (omega > 0.5_dp*sqrt(huge(1.0_dp))) return
        end if
        if (size(a) == 0 .or. size(b) == 0) return
        if (.not. finite_complex(a) .or. .not. finite_complex(b)) return
        if (all(abs(a) == 0.0_dp)) return

        ! Discard leading zero coefficients before choosing the degree.
        ib = 1
        do while (ib < size(b))
            if (abs(b(ib)) /= 0.0_dp) exit
            ib = ib + 1
        end do
        ia = 1
        do while (ia < size(a))
            if (abs(a(ia)) /= 0.0_dp) exit
            ia = ia + 1
        end do
        bc = b(ib:)
        ac = a(ia:)
        bdeg = size(bc)-1
        adeg = size(ac)-1
        degree = max(bdeg, adeg)
        if (degree > max_polynomial_order) then
            out%status = frequency_transform_too_expensive
            return
        end if
        select case (which)
        case (1)  ! t = s/wo
            num = [cmplx(0.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.0_dp,dp)]
            den = [cmplx(omega,0.0_dp,dp)]
        case (2)  ! t = wo/s
            num = [cmplx(omega,0.0_dp,dp)]
            den = [cmplx(0.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.0_dp,dp)]
        case (3)  ! t = (s**2 + wo**2)/(bw*s)
            num = [cmplx(omega**2,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
                cmplx(1.0_dp,0.0_dp,dp)]
            den = [cmplx(0.0_dp,0.0_dp,dp),cmplx(bandwidth,0.0_dp,dp)]
        case (4)  ! t = bw*s/(s**2 + wo**2)
            num = [cmplx(0.0_dp,0.0_dp,dp),cmplx(bandwidth,0.0_dp,dp)]
            den = [cmplx(omega**2,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
                cmplx(1.0_dp,0.0_dp,dp)]
        case default
            return
        end select
        pb = rational_polynomial(bc, degree, num, den)
        pa = rational_polynomial(ac, degree, num, den)
        if (.not. finite_complex(pb) .or. .not. finite_complex(pa)) return

        ! Remove structural leading zeros and make the denominator monic.
        ia = size(pa)
        do while (ia > 1)
            if (abs(pa(ia)) /= 0.0_dp) exit
            ia = ia - 1
        end do
        if (abs(pa(ia)) == 0.0_dp) return
        ib = size(pb)
        do while (ib > 1)
            if (abs(pb(ib)) /= 0.0_dp) exit
            ib = ib - 1
        end do
        scale = pa(ia)
        out%a = pa(ia:1:-1)/scale
        out%b = pb(ib:1:-1)/scale
        if (.not. finite_complex(out%a) .or. .not. finite_complex(out%b)) then
            deallocate(out%a, out%b)
            return
        end if
        out%status = frequency_transform_success
    end function transform_ba

    function lp2lp_real(b, a, wo) result(out)
        real(dp), intent(in) :: b(:) !! real lowpass numerator, descending powers
        real(dp), intent(in) :: a(:) !! real lowpass denominator, descending powers
        real(dp), intent(in), optional :: wo !! new lowpass cutoff, rad/s
        type(analog_tf_transform_result) :: out
        out = transform_ba(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),1,wo)
    end function lp2lp_real

    function lp2hp_real(b, a, wo) result(out)
        real(dp), intent(in) :: b(:) !! real lowpass numerator, descending powers
        real(dp), intent(in) :: a(:) !! real lowpass denominator, descending powers
        real(dp), intent(in), optional :: wo !! highpass corner, rad/s
        type(analog_tf_transform_result) :: out
        out = transform_ba(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),2,wo)
    end function lp2hp_real

    function lp2bp_real(b, a, wo, bw) result(out)
        real(dp), intent(in) :: b(:) !! real lowpass numerator, descending powers
        real(dp), intent(in) :: a(:) !! real lowpass denominator, descending powers
        real(dp), intent(in), optional :: wo !! bandpass center, rad/s
        real(dp), intent(in), optional :: bw !! bandpass bandwidth, rad/s
        type(analog_tf_transform_result) :: out
        out = transform_ba(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),3,wo,bw)
    end function lp2bp_real

    function lp2bs_real(b, a, wo, bw) result(out)
        real(dp), intent(in) :: b(:) !! real lowpass numerator, descending powers
        real(dp), intent(in) :: a(:) !! real lowpass denominator, descending powers
        real(dp), intent(in), optional :: wo !! bandstop center, rad/s
        real(dp), intent(in), optional :: bw !! bandstop bandwidth, rad/s
        type(analog_tf_transform_result) :: out
        out = transform_ba(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),4,wo,bw)
    end function lp2bs_real

    function lp2lp_complex(b, a, wo) result(out)
        complex(dp), intent(in) :: b(:) !! complex numerator, descending powers
        complex(dp), intent(in) :: a(:) !! complex denominator, descending powers
        real(dp), intent(in), optional :: wo !! new lowpass cutoff, rad/s
        type(analog_tf_transform_result) :: out
        out = transform_ba(b,a,1,wo)
    end function lp2lp_complex

    function lp2hp_complex(b, a, wo) result(out)
        complex(dp), intent(in) :: b(:) !! complex numerator, descending powers
        complex(dp), intent(in) :: a(:) !! complex denominator, descending powers
        real(dp), intent(in), optional :: wo !! highpass corner, rad/s
        type(analog_tf_transform_result) :: out
        out = transform_ba(b,a,2,wo)
    end function lp2hp_complex

    function lp2bp_complex(b, a, wo, bw) result(out)
        complex(dp), intent(in) :: b(:) !! complex numerator, descending powers
        complex(dp), intent(in) :: a(:) !! complex denominator, descending powers
        real(dp), intent(in), optional :: wo !! bandpass center, rad/s
        real(dp), intent(in), optional :: bw !! bandpass bandwidth, rad/s
        type(analog_tf_transform_result) :: out
        out = transform_ba(b,a,3,wo,bw)
    end function lp2bp_complex

    function lp2bs_complex(b, a, wo, bw) result(out)
        complex(dp), intent(in) :: b(:) !! complex numerator, descending powers
        complex(dp), intent(in) :: a(:) !! complex denominator, descending powers
        real(dp), intent(in), optional :: wo !! bandstop center, rad/s
        real(dp), intent(in), optional :: bw !! bandstop bandwidth, rad/s
        type(analog_tf_transform_result) :: out
        out = transform_ba(b,a,4,wo,bw)
    end function lp2bs_complex

    function transform_zpk(z, p, k, which, wo, bw) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        complex(dp), intent(in) :: k !! prototype gain
        integer, intent(in) :: which !! 1 lowpass, 2 highpass, 3 bandpass, 4 bandstop
        real(dp), intent(in), optional :: wo !! center frequency, rad/s
        real(dp), intent(in), optional :: bw !! bandwidth, rad/s
        type(analog_zpk_transform_result) :: out
        real(dp) :: omega, bandwidth
        complex(dp) :: gain, v, d, numerator, denominator
        integer :: i, nz, np, excess

        omega = 1.0_dp
        bandwidth = 1.0_dp
        if (present(wo)) omega = wo
        if (present(bw)) bandwidth = bw
        if (.not. ieee_is_finite(omega) .or. .not. ieee_is_finite(bandwidth)) return
        if (omega <= 0.0_dp .or. bandwidth <= 0.0_dp) return
        ! Band transformations explicitly form wo**2; reject an unrepresentable
        ! center-frequency square before evaluating it under FPE traps.
        if (which >= 3) then
            if (omega > 0.5_dp*sqrt(huge(1.0_dp))) return
        end if
        if (.not. finite_complex(z) .or. .not. finite_complex(p)) return
        if (.not. finite_complex([k])) return
        nz = size(z)
        np = size(p)
        if (nz > np) return
        if (np > 2*max_polynomial_order) then
            out%status = frequency_transform_too_expensive
            return
        end if
        excess = np-nz
        gain = k

        select case (which)
        case (1)
            out%zeros = omega*z
            out%poles = omega*p
            gain = gain*omega**excess
        case (2)
            if (any(abs(z) == 0.0_dp) .or. any(abs(p) == 0.0_dp)) return
            allocate(out%zeros(np), out%poles(np))
            out%zeros(:nz) = omega/z
            out%zeros(nz+1:) = (0.0_dp,0.0_dp)
            out%poles = omega/p
            numerator = (1.0_dp,0.0_dp)
            denominator = (1.0_dp,0.0_dp)
            do i = 1, nz
                numerator = numerator*(-z(i))
            end do
            do i = 1, np
                denominator = denominator*(-p(i))
            end do
            gain = gain*numerator/denominator
        case (3)
            allocate(out%zeros(2*nz+excess), out%poles(2*np))
            do i = 1, nz
                v = 0.5_dp*bandwidth*z(i)
                d = sqrt(v*v-cmplx(omega**2,0.0_dp,dp))
                out%zeros(i) = v+d
                out%zeros(i+nz) = v-d
            end do
            out%zeros(2*nz+1:) = (0.0_dp,0.0_dp)
            do i = 1, np
                v = 0.5_dp*bandwidth*p(i)
                d = sqrt(v*v-cmplx(omega**2,0.0_dp,dp))
                out%poles(i) = v+d
                out%poles(i+np) = v-d
            end do
            gain = gain*bandwidth**excess
        case (4)
            if (any(abs(z) == 0.0_dp) .or. any(abs(p) == 0.0_dp)) return
            allocate(out%zeros(2*nz+2*excess), out%poles(2*np))
            do i = 1, nz
                v = 0.5_dp*bandwidth/z(i)
                d = sqrt(v*v-cmplx(omega**2,0.0_dp,dp))
                out%zeros(i) = v+d
                out%zeros(i+nz) = v-d
            end do
            do i = 1, excess
                out%zeros(2*nz+i) = cmplx(0.0_dp,omega,dp)
                out%zeros(2*nz+excess+i) = cmplx(0.0_dp,-omega,dp)
            end do
            do i = 1, np
                v = 0.5_dp*bandwidth/p(i)
                d = sqrt(v*v-cmplx(omega**2,0.0_dp,dp))
                out%poles(i) = v+d
                out%poles(i+np) = v-d
            end do
            numerator = (1.0_dp,0.0_dp)
            denominator = (1.0_dp,0.0_dp)
            do i = 1, nz
                numerator = numerator*(-z(i))
            end do
            do i = 1, np
                denominator = denominator*(-p(i))
            end do
            gain = gain*numerator/denominator
        case default
            return
        end select
        if (.not. finite_complex(out%zeros) .or. .not. finite_complex(out%poles) .or. &
            .not. finite_complex([gain])) then
            deallocate(out%zeros,out%poles)
            return
        end if
        out%gain = gain
        out%status = frequency_transform_success
    end function transform_zpk

    function lp2lp_zpk_real_gain(z, p, k, wo) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        real(dp), intent(in) :: k !! prototype real gain
        real(dp), intent(in), optional :: wo !! center frequency
        type(analog_zpk_transform_result) :: out
        out = transform_zpk(z,p,cmplx(k,0.0_dp,dp),1,wo)
    end function lp2lp_zpk_real_gain

    function lp2hp_zpk_real_gain(z, p, k, wo) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        real(dp), intent(in) :: k !! prototype real gain
        real(dp), intent(in), optional :: wo !! corner frequency
        type(analog_zpk_transform_result) :: out
        out = transform_zpk(z,p,cmplx(k,0.0_dp,dp),2,wo)
    end function lp2hp_zpk_real_gain

    function lp2bp_zpk_real_gain(z, p, k, wo, bw) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        real(dp), intent(in) :: k !! prototype real gain
        real(dp), intent(in), optional :: wo !! center frequency
        real(dp), intent(in), optional :: bw !! bandwidth
        type(analog_zpk_transform_result) :: out
        out = transform_zpk(z,p,cmplx(k,0.0_dp,dp),3,wo,bw)
    end function lp2bp_zpk_real_gain

    function lp2bs_zpk_real_gain(z, p, k, wo, bw) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        real(dp), intent(in) :: k !! prototype real gain
        real(dp), intent(in), optional :: wo !! center frequency
        real(dp), intent(in), optional :: bw !! bandwidth
        type(analog_zpk_transform_result) :: out
        out = transform_zpk(z,p,cmplx(k,0.0_dp,dp),4,wo,bw)
    end function lp2bs_zpk_real_gain

    function lp2lp_zpk_complex_gain(z, p, k, wo) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        complex(dp), intent(in) :: k !! prototype complex gain
        real(dp), intent(in), optional :: wo !! center frequency
        type(analog_zpk_transform_result) :: out
        out = transform_zpk(z,p,k,1,wo)
    end function lp2lp_zpk_complex_gain

    function lp2hp_zpk_complex_gain(z, p, k, wo) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        complex(dp), intent(in) :: k !! prototype complex gain
        real(dp), intent(in), optional :: wo !! corner frequency
        type(analog_zpk_transform_result) :: out
        out = transform_zpk(z,p,k,2,wo)
    end function lp2hp_zpk_complex_gain

    function lp2bp_zpk_complex_gain(z, p, k, wo, bw) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        complex(dp), intent(in) :: k !! prototype complex gain
        real(dp), intent(in), optional :: wo !! center frequency
        real(dp), intent(in), optional :: bw !! bandwidth
        type(analog_zpk_transform_result) :: out
        out = transform_zpk(z,p,k,3,wo,bw)
    end function lp2bp_zpk_complex_gain

    function lp2bs_zpk_complex_gain(z, p, k, wo, bw) result(out)
        complex(dp), intent(in) :: z(:) !! prototype finite zeros
        complex(dp), intent(in) :: p(:) !! prototype finite poles
        complex(dp), intent(in) :: k !! prototype complex gain
        real(dp), intent(in), optional :: wo !! center frequency
        real(dp), intent(in), optional :: bw !! bandwidth
        type(analog_zpk_transform_result) :: out
        out = transform_zpk(z,p,k,4,wo,bw)
    end function lp2bs_zpk_complex_gain

end module scifort_signal_frequency_transforms
